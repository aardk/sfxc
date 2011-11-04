#include "data_reader_blocking.h"
#include "vdif_reader.h"

VDIF_reader::
VDIF_reader(boost::shared_ptr<Data_reader> data_reader,
              Data_frame &data, Time ref_time)
  : Input_data_format_reader(data_reader),
    debug_level_(CHECK_PERIODIC_HEADERS),
    sample_rate(0)
{
  // Reference date : All times are relative to midnight on ref_jday
  ref_jday = (int)ref_time.get_mjd();
  DEBUG_MSG("Ref_jday=" << ref_jday);
}

VDIF_reader::~VDIF_reader() {}

bool 
VDIF_reader::open_input_stream(Data_frame &data){
  if(!read_new_block(data)){
    sfxc_abort("Couldn't find valid VDIF header");
  }

  is_open_ = true;
  epoch_jday = current_header.jday_epoch();
  current_time_ = get_current_time();
  data.start_time = current_time_;
  uint start_sec = current_header.sec_from_epoch;
  uint epoch = current_header.ref_epoch;
  std::cout << RANK_OF_NODE << "Start of VDIF data at jday=" << epoch_jday+start_sec/(24*60*60)
            << ", seconds in epoch = " << start_sec << ", epoch="<< epoch 
            << ", t=" <<  current_time_ << "\n";
  return true;
}

void 
VDIF_reader::print_header(){
  std::cout << RANK_OF_NODE << "------------ full header ------------\n";
  std::cout << current_header.sec_from_epoch<<" ; " << (int)current_header.legacy_mode << " ; "
            << (int)current_header.invalid<<"\n";
  std::cout << current_header.dataframe_in_second << " ; " << (int)current_header.ref_epoch
            << " ; " << (int) current_header.unassiged<<"\n";
  std::cout << current_header.dataframe_length << " ; " << (int) current_header.log2_nchan << " ; "
            << (int)current_header.version<<"\n";
  std::cout << (int) current_header.station_id << " ; " << (int)current_header.thread_id << " ; "
            << (int) current_header.bits_per_sample << " ;" << (int) current_header.data_type<<"\n";
  std::cout << current_header.user_data1<<" ; "<<(int)current_header.edv<<"\n";
  std::cout << current_header.user_data2<<" ; "<<current_header.user_data3<<" ; "<< current_header.user_data4<<"\n";
  std::cout << RANK_OF_NODE << "-------------------------------------\n";
}

Time
VDIF_reader::goto_time(Data_frame &data, Time time) {
  while (current_header.thread_id != 0 && time > get_current_time())
    read_new_block(data);

  return get_current_time();
}

Time VDIF_reader::get_current_time(){
  Time time;

  if (is_open_) {
    double seconds_since_reference = (double)current_header.sec_from_epoch-(ref_jday-epoch_jday)*24*60*60;
    double subsec = 0;
    if(sample_rate>0){
      int samples_per_frame = 8*current_header.data_size()/(current_header.bits_per_sample+1);
      subsec = (double)current_header.dataframe_in_second*samples_per_frame/sample_rate;
    }
    time.set_time(ref_jday, seconds_since_reference + subsec);
  }

  return time;
}

bool VDIF_reader::read_new_block(Data_frame &data) {
  std::vector<value_type> &buffer = data.buffer->data;
  int byte_read = Data_reader_blocking::get_bytes_s( data_reader_.get(), 16, (char *)&current_header );
  if(current_header.legacy_mode==0){
    char *header = (char *)&current_header;
    Data_reader_blocking::get_bytes_s(data_reader_.get(), 16, (char *)&header[16]);
  }
  //  print_header();
  int data_size = current_header.data_size();
  if (buffer.size() != data_size) {
    buffer.resize(data_size);
  }

  if(current_header.invalid > 0){
    std::cout << RANK_OF_NODE << " : invalid frame!!!\n";
    data.invalid.resize(1);
    data.invalid[0].invalid_begin = 0;
    data.invalid[0].nr_invalid = current_header.data_size()*current_header.invalid;
  }
  data.channel = current_header.thread_id; // NB: we assume one-to-one mapping of thread id to channel

  byte_read = Data_reader_blocking::get_bytes_s( data_reader_.get(), data_size, (char *)&buffer[0]);

  if (data_reader_->eof()) return false;

  data.start_time = get_current_time();
  return true;
}

bool VDIF_reader::eof() {
  return data_reader_->eof();
}

int32_t VDIF_reader::Header::jday_epoch() const {
  int year = 2000+ref_epoch/2;
  int month = 1+6*(ref_epoch&1);
  return mjd(1,month,year);
}

Time VDIF_reader::time_between_headers() {
  int samples_per_byte = 8/(current_header.bits_per_sample+1);
  Time time_between_headers_(current_header.data_size() * samples_per_byte / (sample_rate / 1000000));
  SFXC_ASSERT(time_between_headers_.get_time_usec() > 0);
  return time_between_headers_;
}

void VDIF_reader::set_parameters(const Input_node_parameters &param) {
  sample_rate = param.sample_rate();
  SFXC_ASSERT(((int)sample_rate % 1000000) == 0);
}
