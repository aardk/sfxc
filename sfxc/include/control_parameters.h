#ifndef CONTROL_PARAMETERS_H_
#define CONTROL_PARAMETERS_H_

#include <json/json.h>
#include <vex/Vex++.h>
#include <list>
#include <algorithm>
#include <map>
#include <math.h>
#include "utils.h"
#include "correlator_time.h"

typedef std::pair<std::string, std::string> stream_key;

/** Information about the mark5 tracks needed by the input node. **/
class Input_node_parameters {
public:
  Input_node_parameters()
      : track_bit_rate(0), data_modulation(0) {}

  class Channel_parameters {
  public:
    Channel_parameters() : bits_per_sample(-1) {}


    bool operator==(const Channel_parameters &other) const;

    int bits_per_sample;                   ///< Number of bits to encode one sample
    std::vector<int32_t> tracks;           ///< A list of bit positions for all tracks

    char sideband;
    char polarisation;
    int32_t frequency_number;
    int32_t extra_delay_in_samples;
  };

  typedef std::vector<Channel_parameters>           Channel_list;
  typedef Channel_list::iterator                    Channel_iterator;
  typedef Channel_list::const_iterator              Channel_const_iterator;

  int bits_per_sample() const;
  int subsamples_per_sample() const;
  uint64_t sample_rate() const;
  bool operator==(const Input_node_parameters &other) const;

  /// List of the tracks that are combined to frequency channels
  Channel_list channels;
  /// Number of tracks / bitstreams in data file
  int32_t n_tracks;
  // data rate of one subband
  uint64_t track_bit_rate; // in Ms/s
  // Frame size
  int32_t frame_size;
  /// Time offset for the data reader, used e.g. to compensate for formatter errors
  Time offset;
  /// Indicates if data modulation is used (p.6 of Mark4 memo 230A, Whitney 2005)
  int32_t data_modulation;
  // Phasecal integration time
  Time phasecal_integr_time;
  // Because of dispersion, consecutive time slices overlap
  Time overlap_time;
  // Abort the correlation if the input stream contains no valid data
  bool exit_on_empty_datastream;
};

std::ostream &operator<<(std::ostream &out, const Input_node_parameters &param);


class Pulsar_parameters {
public:
  Pulsar_parameters(std::ostream& log_writer_);
  /// All parameters from the polyco file
  struct Polyco_params{
    char name[11];            // Pulsar name
    char date[10];            // Date 
    float utc;                // UTC 
    double tmid;              // Reference time [mjd]
    double DM;                // Dispersion measure
    float doppler;            // Doppler shift due to earth motion (10^-4)
    float residual;           // Log10 of fit rms residual in periods
    double ref_phase;         // Reference phase (RPHASE)
    double ref_freq;          // Reference rotation frequency (F0)
    char site[6];             // Observatory name
    int32_t data_span;        // Data span [minutes]
    int32_t n_coef;           // The number of coefficients in the polynomial
    double obs_freq;          // Observing frequency
    float bin_phase[2];        // Binary phase
    std::vector<double> coef; // The polynomial coefficients
  };

  struct Pulsar{
    char name[11];
    int32_t nbins;
    struct Interval{double start; double stop;} interval;
    std::vector<Polyco_params> polyco_params;
  };

  bool parse_polyco(std::vector<Polyco_params> &param, std::string filename);

  // maps pulsar name to vector of pulsar pa
  std::map<std::string, Pulsar> pulsars;
private:
  std::ostream& log_writer;
};

class Mask_parameters {
 public:
   Mask_parameters() : normalize(false) {}

  bool normalize;
  std::vector<double> mask;
  std::vector<double> window;
};

/** Information about the correlation neede by the correlator node. **/
class Correlation_parameters {
public:
  Correlation_parameters()
    : number_channels(0), fft_size_delaycor(0), fft_size_correlation(0),
    fft_size_dedispersion(0), integration_nr(-1), slice_nr(-1), sample_rate(0),
    channel_freq(0), bandwidth(0), sideband('n'), frequency_nr(-1),
    polarisation('n'), multi_phase_center(false), pulsar_binning(false),
    window(SFXC_WINDOW_RECT) {}

  bool operator==(const Correlation_parameters& other) const;

  class Station_parameters {
  public:
    bool
    operator==(const Correlation_parameters::Station_parameters& other) const;

    int32_t station_number; // the number of the station
    // according to the vex file
    // sorted alphabathically
    int32_t station_stream; // input stream (from multiple_data_readers)
    uint64_t sample_rate;
    int64_t channel_freq;
    uint64_t bandwidth;
    char sideband;
    char polarisation;
    int32_t bits_per_sample;
    double LO_offset; // LO offset in Hz
    double extra_delay;
    int tsys_freq;
  };

  typedef std::vector<Station_parameters> Station_list;
  typedef Station_list::iterator          Station_iterator;

  // Data members
  Time experiment_start;    // Start time of the experiment
  int64_t slice_size;       // Number of samples in slice
  Time integration_start;   // Start of the integration
  Time integration_time;    // Length of an integration
  Time slice_start;	    // Start of the integration slice
  Time slice_time;          // Length of an integration slice
  Time sub_integration_time;// Length of a subintegration
  Time stream_start;        // Start of the data (for dedispersion)
  int32_t number_channels;  // number of frequency channels
  int32_t fft_size_delaycor;    // Number of samples per FFT in the delay correction
  int32_t fft_size_dedispersion; // Number of samples per FFT in the coherent dedispersion
  int32_t fft_size_correlation; // Number of samples per FFT in the (cross-)correlation
  int32_t integration_nr;   // number of the integration
  int32_t slice_nr;         // Number of the output slice
  // between one integration slice and the next
  // in case of subsecond integrations
  uint64_t sample_rate;     // #Samples per second
  int64_t channel_freq;     // Center frequency of the band in Hz
  uint64_t bandwidth;       // Bandwidth of the channel in Hz
  char    sideband;         // U or L
  int32_t frequency_nr;     // Canonical frequency number
  char    polarisation;     // L or R

  bool    cross_polarize;   // do the cross polarisations
  int32_t reference_station;// use a reference station

  Station_list station_streams; // input streams used
  int window;                   // Windowing function to be used
  char source[11];              // name of the source under observation
  int32_t n_phase_centers;   // The number of phase centers in the current scan
  int32_t multi_phase_center;
  int32_t pulsar_binning;
  Pulsar_parameters *pulsar_parameters;
  Mask_parameters *mask_parameters;
};


std::ostream &operator<<(std::ostream &out, const Correlation_parameters &param);

/** Class containing all control variables needed for the experiment **/
class Control_parameters {
public:

  Control_parameters();
  Control_parameters(const char *ctrl_file, const char *vex_file,
                     std::ostream& log_writer);

  bool initialise(const char *ctrl_filename,
                  const char *vex_filename,
                  std::ostream& log_writer);

  bool check(std::ostream &log_writer) const;

  bool get_pulsar_parameters(Pulsar_parameters &pars) const;
  bool get_mask_parameters(Mask_parameters &pars) const;

  /****************************************************/
  /* Get functions from the correlation control file: */
  /****************************************************/
  Time get_start_time() const;
  Time get_stop_time() const;
  std::vector<std::string> data_sources(const std::string &station) const;
  std::vector<std::string> data_sources(const std::string &station,
					const std::string &datastream) const;
  std::string get_output_file() const;
  std::string get_phasecal_file() const;
  std::string get_tsys_file() const;

  std::string station(int i) const;
  size_t number_stations() const;
  int station_number(const std::string &station_name) const;
  size_t number_inputs() const;

  Time integration_time() const; // Integration time in microseconds
  Time sub_integration_time() const;
  Time phasecal_integration_time() const;
  int slices_per_integration() const;
  int number_channels() const;
  int fft_size_delaycor() const;
  int fft_size_correlation() const;
  int window_function() const;
  int job_nr() const;
  int subjob_nr() const;
  int output_buffer_size() const;

  std::string sideband(int i) const;
  std::string reference_station() const;
  int reference_station_number() const;
  std::string setup_station() const;

  bool phased_array() const;
  bool pulsar_binning() const;
  bool multi_phase_center() const;
  double LO_offset(const std::string &station, int integration_nr) const;
  double extra_delay(const std::string &channel_name,
		     const std::string &station_name,
		     const std::string &mode_name) const;
  int extra_delay_in_samples(const std::string &channel_name,
			     const std::string &station_name,
			     const std::string &mode_name) const;
  int tsys_freq(const std::string &station) const;
  bool exit_on_empty_datastream() const;
  
  Time reader_offset(const std::string &s) const{
    return reader_offsets.find(s)->second;
  };
  
  void set_reader_offset(const std::string &s, const Time t);

  std::string get_delay_table_name(const std::string &station_name) const;
  void generate_delay_table(const std::string &station_name,
                            const std::string &filename) const;
  std::string channel(int i) const;

  int message_level() const;

  /****************************************************/
  /* Get functions from the vex file:                 */
  /****************************************************/
  int bits_per_sample(const std::string& mode, const std::string& station) const;
  uint64_t sample_rate(const std::string& mode, const std::string& station) const;
  uint64_t bandwidth(const std::string& mode, const std::string& station, const std::string& channel) const;
  int64_t channel_freq(const std::string& mode, const std::string& station, const std::string& channel) const;
  std::string datastream(const std::string& mode, const std::string& station, const std::string& channel) const;
  std::vector<std::string> datastreams(const std::string& station) const;

  std::string scan(int i) const;
  int scan(const Time &time) const;
  std::string scan_source(const std::string &scan) const;

  size_t number_scans() const;

  bool station_in_scan(const std::string& scan, const std::string &station) const;
  size_t number_stations_in_scan(const std::string& scan) const;

  Time stop_time(const std::string& scan, const std::string &station) const;

  // Takes cross polarisation into account
  int number_correlation_cores_per_timeslice(const std::string &mode) const;

  // Return the Frequency channels from the VEX file, filtered by the ctrl file
  size_t number_frequency_channels() const;
  std::string frequency_channel(size_t channel_nr, const std::string& mode_name, const std::string &station_name) const;
  int frequency_number(size_t channel_nr, const std::string& mode_name) const;

  bool cross_polarize() const;
  int cross_channel(int channel_nr,
                    const std::string &mode) const;
  int cross_channel(const std::string &channel_nr,
                    const std::string &mode) const;

  char polarisation(const std::string &channel_name,
                    const std::string &station_name,
                    const std::string &mode) const;

  std::string frequency(const std::string &channel_name,
                        const std::string &station_name,
                        const std::string &mode) const;

  char sideband(const std::string &channel_name,
                const std::string &station_name,
                const std::string &mode) const;

  // Returns the number of correlation ffts to be processed for one integration slice
  static int nr_correlation_ffts_per_integration(const Time &integration_time,
				           uint64_t sample_rate, int fft_size) {
    Time time_one_fft(fft_size / (sample_rate / 1000000.));
    return floor(integration_time / time_one_fft);
  }

  // Returns the number of de-dispersion ffts to be processed for one integration slice, 
  // taking into account that integrations are half overlapped
  static int nr_dedisp_ffts_per_integration(const Time &integration_time,
                                         int sample_rate,
                                         int fft_size_dedispersion,
                                         int fft_size_correlation) {
    Time time_correlation(fft_size_correlation / (sample_rate / 1000000.));
    Time time_dedispersion(fft_size_dedispersion / (sample_rate / 1000000.));
    int nr_corr_fft = nr_correlation_ffts_per_integration(integration_time,
                        sample_rate, fft_size_correlation);
    int nr_fft = ceil((time_correlation*nr_corr_fft + time_dedispersion) /
                      time_dedispersion);
    return nr_fft;
  }

  int polarisation_type_for_global_output_header(const std::string &mode) const;

  /****************************************************/
  /* Extract structs for the correlation:             */
  /****************************************************/

  // Return the track parameters needed by the input node
  Input_node_parameters
  get_input_node_parameters(const std::string &mode_name,
                            const std::string &station_name,
			    const std::string &datastream_name) const;

  // Return the correlation parameters needed by a correlator node
  Correlation_parameters
  get_correlation_parameters(const std::string &scan_name,
                             size_t channel_nr,
                             int integration_nr,
                             const std::map<stream_key, int> &correlator_node_station_to_input) const;
  std::string rack_type(const std::string &station) const;
  std::string recorder_type(const std::string &station) const;
  std::string data_format(const std::string &station, const std::string &mode) const;

  const Vex &get_vex() const;
  std::string get_exper_name() const;

  /// The start time of the first scan in the vex file
  Time start_time;

private:
  std::string create_path(const std::string &path) const;
private:

  // Gets the track parameters for mark5a data
  // Output is in input_parameters
  void get_mark5a_tracks(const std::string &mode,
                         const std::string &station,
                         Input_node_parameters &input_parameters) const;
  // Get the bit positions for all tracks in the vex file
  std::vector<int> get_track_bit_position(const std::string &mode, const std::string &station) const;
  int n_mark5a_tracks(const std::string &mode, const std::string &station) const;

  // Gets the track parameters for mark5b data
  // Output is in input_parameters
  void get_mark5b_tracks(const std::string &mode,
                         const std::string &station,
                         Input_node_parameters &input_parameters) const;
  int n_mark5b_bitstreams(const std::string &mode, const std::string &station) const;
  void get_mark5b_standard_mapping(const std::string &mode,
                                   const std::string &station,
                                   Input_node_parameters &input_parameters) const;

  // Gets the track parameters for VDIF data
  // Output is in input_parameters
  void get_vdif_tracks(const std::string &mode,
		       const std::string &station,
		       const std::string &datastream,
		       Input_node_parameters &input_parameters) const;
  void get_vdif_datastreams(const std::string &mode,
			    const std::string &station,
			    const std::string &datastream,
			    Input_node_parameters &input_parameters) const;
  void get_vdif_threads(const std::string &mode,
			const std::string &station,
			Input_node_parameters &input_parameters) const;

  std::string ctrl_filename;
  std::string vex_filename;
  std::map<std::string, Time> reader_offsets; // Contains the formatter clock offsets for all input nodes

  Json::Value ctrl;        // Correlator control file
  Vex         vex;         // Vex file
  bool        initialised; // The control parameters are initialised

  mutable std::map<std::string, int> station_map;

  bool check_data_source(std::ostream &log_writer, const Json::Value &) const;
};

#endif /*CONTROL_PARAMETERS_H_*/
