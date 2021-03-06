/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *            Aard Keimpema <Keimpema@JIVE.nl>, 2008
 *
 * $Id: correlator_node.cc 2215 2019-09-04 10:10:54Z kettenis $
 *
 */

#include "correlator_node_tasklet.h"
#include "correlation_core_phased.h"
#include "data_reader_buffer.h"
#include "data_writer.h"
#include "utils.h"
#include "output_header.h"
#include "delay_correction.h"

Correlator_node_tasklet::Correlator_node_tasklet(int nr_corr_node, bool pulsar_binning_, bool phased_array_) :
    status(STOPPED),
    isinitialized_(false),
    nr_corr_node(nr_corr_node),
    pulsar_binning(pulsar_binning_),
    phased_array(phased_array_),
    has_requested(false) {
  if (phased_array){
    correlation_core_normal = new Correlation_core_phased();
    correlation_core = correlation_core_normal;
  }else{
    correlation_core_normal = new Correlation_core();
    if(pulsar_binning){
      correlation_core_pulsar = new Correlation_core_pulsar();
      correlation_core = correlation_core_pulsar;
    }else
      correlation_core = correlation_core_normal;
  }

#ifdef RUNTIME_STATISTIC

  std::stringstream inputid;
  std::stringstream compid;
  std::stringstream monid;
  std::stringstream tt;

  inputid << "correlationnode" << RANK_OF_NODE;

  compid << inputid.str() << "_dotask";
  monid << compid.str() << "_monitor_state";
  dotask_state_.init(monid.str());
  dotask_state_.add_property(inputid.str(), "is_a", "correlationnode");
  dotask_state_.add_property(inputid.str(), "has", compid.str() );
  dotask_state_.add_property(compid.str(), "is_a", "correlationnode_dotaskloop");
  dotask_state_.add_property(compid.str(), "has", monid.str() );
  tt.str(monid.str());

  compid.str("");
  monid.str("");
  compid << inputid.str() << "_reader";
  monid << compid.str() << "_monitor_state";
  reader_state_.init(monid.str());
  reader_state_.add_property(inputid.str(), "is_a", "correlationnode");
  reader_state_.add_property(inputid.str(), "has", compid.str() );
  reader_state_.add_property(compid.str(), "is_a", "reader_component");
  reader_state_.add_property(compid.str(), "has", monid.str() );
  dotask_state_.add_property(tt.str(), "contains", monid.str() );

  compid.str("");
  monid.str("");
  compid << inputid.str() << "_delaycorrection";
  monid << compid.str() << "_monitor_state";
  delaycorrection_state_.init(monid.str());
  delaycorrection_state_.add_property(inputid.str(), "is_a", "correlationnode");
  delaycorrection_state_.add_property(inputid.str(), "has", compid.str() );
  delaycorrection_state_.add_property(compid.str(), "is_a", "delaycorrection_component");
  delaycorrection_state_.add_property(compid.str(), "has", monid.str() );
  dotask_state_.add_property(tt.str(), "contains", monid.str() );

  compid.str("");
  monid.str("");
  compid << inputid.str() << "_integration";
  monid << compid.str() << "_monitor_state";
  correlation_state_.init(monid.str());
  correlation_state_.add_property(inputid.str(), "is_a", "correlationnode");
  correlation_state_.add_property(inputid.str(), "has", compid.str() );
  correlation_state_.add_property(compid.str(), "is_a", "correlation_component");
  correlation_state_.add_property(compid.str(), "has", monid.str() );
  correlation_state_.add_property(tt.str(), "contains", monid.str() );
#endif //RUNTIME_STATISTIC
}

Correlator_node_tasklet::~Correlator_node_tasklet() {
#if PRINT_TIMER
  PROGRESS_MSG("Time bit_sample_reader:  " << bit_sample_reader_timer_.measured_time());
  PROGRESS_MSG("Time bits2float:  " << bits_to_float_timer_.measured_time());
  PROGRESS_MSG("Time delay:       " << delay_timer_.measured_time());
  PROGRESS_MSG("Time correlation: " << correlation_timer_.measured_time());
#endif
}

void Correlator_node_tasklet::start_threads() {
  threadpool_.register_thread( reader_thread_.start() );
  threadpool_.register_thread( bit2float_thread_.start() );
}

void Correlator_node_tasklet::stop_threads() {
  reader_thread_.stop();
  bit2float_thread_.stop();

  /// We wait the termination of the threads
  threadpool_.wait_for_all_termination();
}

void Correlator_node_tasklet::do_execute() {
  // Request first slice
  MPI_Send(&nr_corr_node, 1, MPI_INT32,
           RANK_MANAGER_NODE, MPI_TAG_CORRELATION_OF_TIME_SLICE_ENDED,
           MPI_COMM_WORLD);

  has_requested = true;

  try {
    while (isrunning_) {
      switch (status) {
      case STOPPED: {
        // blocking:
        const Correlation_parameters &parameters = integration_slices_queue.front();
        set_parameters(parameters);
        integration_slices_queue.pop();
        break;
      }
      case CORRELATING: {
        correlate();
        if (!has_requested && correlation_core->almost_finished()) {
          int32_t msg = get_correlate_node_number();
          MPI_Send(&msg, 1, MPI_INT32, RANK_MANAGER_NODE,
                   MPI_TAG_CORRELATION_OF_TIME_SLICE_ENDED,
                   MPI_COMM_WORLD);

          has_requested = true;
        }
        if (correlation_core->finished()) {
          status = STOPPED;
        }
        break;
      }
    } 
    }
  } catch (QueueClosedException &e) {
    /// Integration slice queue is closed: Stop correlation
  }
  stop_threads();
}

void Correlator_node_tasklet::terminate() {
  isrunning_ = false;
  integration_slices_queue.close();
}

void Correlator_node_tasklet::add_delay_table(Delay_table &table, int sn1, int sn2) {
  SFXC_ASSERT(sn1<=sn2);
  SFXC_ASSERT((size_t)sn1 < delay_modules.size());
  SFXC_ASSERT((size_t)sn2 < delay_modules.size());
  SFXC_ASSERT(delay_modules[sn1] != Delay_correction_ptr());
  SFXC_ASSERT(delay_modules[sn2] != Delay_correction_ptr());
  if(delay_tables.size() <= sn1)
    delay_tables.resize(sn1+1);
  if(delay_index.size() <= sn2)
    delay_index.resize(sn2+1, -1);

  delay_tables[sn1].add_scans(table);
  delay_index[sn1] = sn1;
  delay_index[sn2] = sn1;
}

void Correlator_node_tasklet::add_uvw_table(Uvw_model &table, int sn) {
  if(uvw_tables.size() <= sn){
    uvw_tables.resize(sn+1);
  }
  uvw_tables[sn].add_scans(table);
}

void Correlator_node_tasklet::hook_added_data_reader(size_t stream_nr, Data_reader_ptr data_reader) {
  // create the bit sample reader tasklet
  if (reader_thread_.bit_sample_readers().size() <= stream_nr) {
    reader_thread_.bit_sample_readers().resize(stream_nr+1, Bit_sample_reader_ptr());
  }
  reader_thread_.bit_sample_readers()[stream_nr] =
       Bit_sample_reader_ptr(new Correlator_node_data_reader_tasklet());
  reader_thread_.bit_sample_readers()[stream_nr]->connect_to(stream_nr, data_reader);

  // connect reader to data stream worker

  bit_statistics_ptr statistics = bit_statistics_ptr(new bit_statistics());
  bit2float_thread_.connect_to(stream_nr, statistics,
                               reader_thread_.bit_sample_readers()[stream_nr]->get_output_buffer());

  { // create the delay modules
    if (delay_modules.size() <= stream_nr) {
      delay_modules.resize(stream_nr+1, shared_ptr<Delay_correction>());
    }
    delay_modules[stream_nr] = Delay_correction_ptr(new Delay_correction(stream_nr));
    // Connect the delay_correction to the bits2float_converter
    delay_modules[stream_nr]->connect_to(bit2float_thread_.get_output_buffer(stream_nr));
  }


  // Connect the correlation_core to delay_correction
  correlation_core_normal->connect_to(stream_nr, statistics,
                                      delay_modules[stream_nr]->get_output_buffer());
  correlation_core_normal->connect_to(stream_nr, bit2float_thread_.get_invalid(stream_nr));
  if(pulsar_binning){
    correlation_core_pulsar->connect_to(stream_nr, statistics,
                                        delay_modules[stream_nr]->get_output_buffer());
    correlation_core_pulsar->connect_to(stream_nr, bit2float_thread_.get_invalid(stream_nr));
  }
}

void Correlator_node_tasklet::hook_added_data_writer(size_t i, Data_writer_ptr data_writer) {
  SFXC_ASSERT(i == 0);

  correlation_core_normal->set_data_writer(data_writer);
  if(pulsar_binning)
    correlation_core_pulsar->set_data_writer(data_writer);
}

int Correlator_node_tasklet::get_correlate_node_number() {
  return nr_corr_node;
}

void Correlator_node_tasklet::correlate() {
  RT_STAT( dotask_state_.begin_measure() );
  bool done_work=false; 
  delay_timer_.resume();
  for (size_t i=0; i<delay_modules.size(); i++) {
    if (delay_modules[i] != Delay_correction_ptr()) {
      if (delay_modules[i]->has_work()) {
        RT_STAT( delaycorrection_state_.begin_measure() );
        delay_modules[i]->do_task();
        RT_STAT( delaycorrection_state_.end_measure(1) );
        done_work=true;
      }
    }
  }
  delay_timer_.stop();

  correlation_timer_.resume();
  if (correlation_core->has_work()) {
    RT_STAT( correlation_state_.begin_measure() );

    correlation_core->do_task();
    done_work=true;

    RT_STAT( correlation_state_.end_measure(1) );
  }

  correlation_timer_.stop();

  RT_STAT( dotask_state_.end_measure(1) );

  if (!done_work)
    usleep(1000);

}

void
Correlator_node_tasklet::add_new_slice(const Correlation_parameters &parameters) {

  integration_slices_queue.push(parameters);

  /// We add the new timeslice to the readers.
  reader_thread_.add_time_slice_to_read(parameters);
}

void
Correlator_node_tasklet::add_source_list(const std::map<std::string, int> &sources) {
  correlation_core_normal->add_source_list(sources);
}

void
Correlator_node_tasklet::set_parameters(const Correlation_parameters &parameters) {
  SFXC_ASSERT(status == STOPPED);

  if ( !isinitialized_ ) {
    ///DEBUG_MSG("START THE THREADS !");
    isinitialized_ = true;
    start_threads();
  }


  // Get delay and UVW tables
  // NB: The total number of streams in the job is not nessecarily the same 
  // as the number of streams in the scan
  Time tmid = parameters.integration_start + parameters.integration_time / 2;
  std::vector<Delay_table_akima> akima_tables(delay_index.size());
  std::vector<std::vector<double> > uvw(uvw_tables.size());
  for(int i=0;i<parameters.station_streams.size();i++){
    int stream = parameters.station_streams[i].station_stream;
    int index = delay_index[stream];
    SFXC_ASSERT(index != -1);
    akima_tables[stream] = 
       delay_tables[index].create_akima_spline(parameters.integration_start,
                                                parameters.integration_time);
    if(stream < uvw.size()){
      uvw[stream].resize(parameters.n_phase_centers*3);
      for(int j=0;j<parameters.n_phase_centers;j++){
        double *out = &uvw[stream][3*j];
        uvw_tables[stream].get_uvw(j, tmid, &out[0], &out[1], &out[2]);
      }
    }
  }
  int nBins=1;
  if(pulsar_binning){
    Pulsar_parameters *pulsar_parameters = parameters.pulsar_parameters;
    std::map<std::string, Pulsar_parameters::Pulsar>::iterator cur_pulsar_it =
                           pulsar_parameters->pulsars.find(std::string(&parameters.source[0]));
    if(cur_pulsar_it == pulsar_parameters->pulsars.end()){
      // Current source is not a pulsar
      nBins = 1;
      correlation_core = correlation_core_normal;
      correlation_core->set_parameters(parameters, akima_tables, uvw, get_correlate_node_number());
    }else{
      Pulsar_parameters::Pulsar &pulsar = cur_pulsar_it->second;
      nBins = pulsar.nbins + 1; // One extra for off-pulse data 
      correlation_core = correlation_core_pulsar;
      correlation_core_pulsar->set_parameters(parameters, pulsar, akima_tables, uvw, get_correlate_node_number());
    }
  }else{
    nBins = parameters.n_phase_centers;
    correlation_core->set_parameters(parameters, akima_tables, uvw, get_correlate_node_number());
  }

  for (size_t i=0; i<delay_modules.size(); i++) {
    if (delay_modules[i] != Delay_correction_ptr()) {
      delay_modules[i]->set_parameters(parameters, akima_tables[i]);
    }
  }
  bit2float_thread_.set_parameters(parameters, akima_tables);

  has_requested=false;
  status = CORRELATING;

  // set the output stream
  int nstreams = parameters.station_streams.size();
  std::set<int> stations_set;
  for (int i = 0; i < nstreams; i++) {
    int station = parameters.station_streams[i].station_number;
    stations_set.insert(station);
  }
  int nstations = stations_set.size();
  int nBaselines = correlation_core->number_of_baselines();
  int size_of_one_baseline = sizeof(std::complex<float>) * (parameters.number_channels + 1);

  int size_uvw = nstations*sizeof(Output_uvw_coordinates);
  // when the cross_polarize flag is set then the correlator node receives 2 polarizations
  int size_stats = nstreams * sizeof(Output_header_bitstatistics);

  int band = parameters.frequency_nr << 2;
  if (parameters.sideband == 'U')
    band |= 2;
  if (parameters.polarisation == 'L')
    band |= 1;

  bool accum = true;
  if (parameters.slice_start + parameters.slice_time >=
      parameters.integration_start + parameters.integration_time)
    accum = false;

  int slice_size;
  slice_size = sizeof(int32_t) + sizeof(Output_header_timeslice) + size_uvw + size_stats +
               nBaselines * ( size_of_one_baseline + sizeof(Output_header_baseline));
  SFXC_ASSERT(nBins >= 1);

  output_node_set_timeslice(parameters.slice_nr, get_correlate_node_number(),
                            band, accum, slice_size, nBins);
}

void
Correlator_node_tasklet::
output_node_set_timeslice(int slice_nr, int stream_nr, int band, int accum,
                          int bytes, int bins) {
  correlation_core->data_writer()->set_size_dataslice(bins * bytes);
  int32_t msg_output_node[] = {stream_nr, slice_nr, band, accum, bytes, bins};
  MPI_Send(&msg_output_node, 6, MPI_INT32,
           RANK_OUTPUT_NODE,
           MPI_TAG_OUTPUT_STREAM_SLICE_SET_ORDER,
           MPI_COMM_WORLD);
}

void Correlator_node_tasklet::get_state(std::ostream &out) {
  out << "\t\"Correlator_node_tasklet\": {\n"
      << "\t\t\"nr_corr_node\": " << nr_corr_node << ",\n"
      << "\t\t\"pulsar_binning\": " << std::boolalpha << pulsar_binning << ",\n"
      << "\t\t\"phased_array\": " << std::boolalpha << phased_array << ",\n"
      << "\t\t\"has_requested\": " << std::boolalpha << has_requested << ",\n"
      << "\t\t\"isrunning\": " << std::boolalpha << isrunning_ << ",\n"
      << "\t\t\"n_integration_slices\": " << integration_slices_queue.size() << ",\n"
      << "\t\t\"state\": ";
  switch (status) {
      case STOPPED: {
        out << "\"STOPPED\"\n";
        break;
      } case CORRELATING: {
        out << "\"CORRELATING\"\n";
        break;
      } default: {
        out << "\"UNKNOWN STATE\"\n";
      }
  }
  out << "\t},\n";
  reader_thread_.get_state(out);
  bit2float_thread_.get_state(out);
  out << "\t\"Delay_correction\": [\n";
  for (size_t i=0; i<delay_modules.size(); i++) {
    if (delay_modules[i] != Delay_correction_ptr()) {
      delay_modules[i]->get_state(out);
    } else {
      out << "\t\t{}";
    }
    if (i < delay_modules.size() - 1)
      out << ",\n";
    else
      out << "\n";
  }
  out << "\t],\n";
  correlation_core->get_state(out);
}
