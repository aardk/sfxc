#ifndef CORRELATION_CORE_VOLTAGES_H_
#define CORRELATION_CORE_VOLTAGES_H_

#include "correlation_core.h"
#include "bandpass.h"
typedef Correlator_node_types::Channel_queue         Channel_queue;
typedef Correlator_node_types::Channel_queue_ptr     Channel_queue_ptr;

class Correlation_core_voltages : public Correlation_core {
public:
  Correlation_core_voltages();
  virtual ~Correlation_core_voltages();
  virtual void do_task();
  virtual bool has_work();
  virtual void set_parameters(const Correlation_parameters &parameters,
                              std::vector<std::vector<double> > &uvw,
                              int node_nr);
  using Correlation_core::connect_to;
  virtual void connect_to(size_t stream, Channel_queue_ptr buffer);
protected:
  void integration_initialise();
  void integration_step(int ostart, int istart, int index);
  void create_baselines(const Correlation_parameters &parameters);
  void integration_write_subints(std::vector< std::vector<float> > &integration_buffer);
  std::vector<Correlator_node_types::Channel_queue_ptr> input_buffers;
  std::vector< FLOAT* >                   input_elements;
  std::vector< std::vector< std::vector<float> > >  station_voltage;
private:
  int ffts_to_skip;
};
#endif /*CORRELATION_CORE_H_*/
