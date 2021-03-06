#include "correlator_node_data_reader_tasklet.h"

Correlator_node_data_reader_tasklet::
Correlator_node_data_reader_tasklet()
  : input_buffer(37100000), bytes_left(0), stream_nr(-1),
    new_stream_available(false),state(IDLE) {
}

Correlator_node_data_reader_tasklet::
~Correlator_node_data_reader_tasklet() {}

/// Set the input
void
Correlator_node_data_reader_tasklet::
connect_to(int stream_nr_, Data_reader_ptr reader_) {
  stream_nr = stream_nr_;
  reader = reader_;
  breader_ = Data_reader_blocking_ptr( new Data_reader_blocking( reader_.get() ) );
}

Correlator_node_data_reader_tasklet::Input_buffer_ptr
Correlator_node_data_reader_tasklet::
get_output_buffer() {
  return &input_buffer;
}

void
Correlator_node_data_reader_tasklet::do_task() {
  uint8_t header;
  std::vector<unsigned char> &data = input_buffer.data;
  size_t dsize = data.size();

  uint64_t read = input_buffer.read;
  uint64_t write = input_buffer.write;

  SFXC_ASSERT(read <= write);

  switch (state) {
  case IDLE:
    if (!new_stream_available)
      return;

    new_stream_available = false;
    state = PROCESSING_STREAM;
    /* FALLTHROUGH */
  case PROCESSING_STREAM:
    breader_->get_bytes(sizeof(header), (char *)&header);
    SFXC_ASSERT(write < (read + dsize - 3));
    data[write++ % dsize] = header;
    switch (header) {
    case HEADER_DATA:
    {
      uint16_t nbytes;
      breader_->get_bytes(sizeof(nbytes), (char *)&nbytes);
      data[write++ % dsize] = nbytes & 0xff;
      data[write++ % dsize] = nbytes >> 8;
      bytes_left = nbytes;
      SFXC_ASSERT(nbytes > 0);
      state = RECEIVE_DATA;
      break;
    }
    case HEADER_DELAY:
    {
      int8_t new_delay;
      breader_->get_bytes(sizeof(new_delay), (char *)&new_delay);
      data[write++ % dsize] = new_delay;
      break;
    }
    case HEADER_INVALID:
    {
      uint16_t n_invalid;
      breader_->get_bytes(sizeof(n_invalid), (char *)&n_invalid);
      data[write++ % dsize] = n_invalid & 0xff;
      data[write++ % dsize] = n_invalid >> 8;
      break;
    }
    case HEADER_ENDSTREAM:
      state = IDLE;
      break;
    default:
      SFXC_ASSERT_MSG(false, "Read invalid header from data stream");
    }
    break;
  case RECEIVE_DATA:
    if (bytes_left > 0) {
      size_t bytes_left_in_buffer = dsize + (read - write);
      size_t to_read = std::min(bytes_left, bytes_left_in_buffer - 1);
      size_t data_read = 0;
      while (data_read < to_read) {
        size_t data_to_read = std::min(to_read - data_read, (size_t)(data.size() - (write % dsize)));
        SFXC_ASSERT(data_to_read > 0);
        size_t nbytes = breader_->get_bytes(data_to_read, (char *)&data[write % dsize]);
        SFXC_ASSERT(nbytes > 0);
        data_read += nbytes;
        write += nbytes;
      }
      bytes_left -= to_read;
    }
    if (bytes_left == 0) {
      state = PROCESSING_STREAM;
    }
    SFXC_ASSERT(bytes_left >= 0);
    break;
  }

  SFXC_ASSERT(read <= write);
  input_buffer.write = write;
}

bool
Correlator_node_data_reader_tasklet::has_work() {
  if (reader == Data_reader_ptr())
    return false;

  if (!reader->can_read())
    return false;

  if(input_buffer.bytes_free() < INPUT_BUFFER_MINIMUM_FREE)
    return false;

  if (state == IDLE)
    return new_stream_available;

  return true;
}

int Correlator_node_data_reader_tasklet::get_fd() {
  return reader->get_fd();
}

bool Correlator_node_data_reader_tasklet::active() {
  if(state!=IDLE)
    return true;
  else
    return new_stream_available;
}

void
Correlator_node_data_reader_tasklet::set_parameters() {
  new_stream_available = true;
}

void Correlator_node_data_reader_tasklet::get_state(std::ostream &out) {
  out << "\t\t{\n"
      << "\t\t\"stream_nr\": " << stream_nr << ",\n"
      << "\t\t\"read\": " << input_buffer.read << ",\n"
      << "\t\t\"write\": " << input_buffer.read << ",\n"
      << "\t\t\"new_stream_available\": " << new_stream_available << ",\n"
      << "\t\t\"state\": ";
  switch (state) {
  case IDLE:
    out << "\"IDLE\"\n";
    break;
  case PROCESSING_STREAM:
    out << "\"PROCESSING_STREAM\"\n";
    break;
  case RECEIVE_DATA:
    out << "\"RECEIVE_DATA\"\n";
    break;
  default:
    out << "\"UNKNOWN_STATE\"\n";
  }
  out << "\t\t}";
}
