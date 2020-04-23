#ifndef CORRELATION_CORE_FILTERBANK_H_
#define CORRELATION_CORE_FILTERBANK_H_

#include "correlation_core.h"
#include "bandpass.h"

class Correlation_core_filterbank : public Correlation_core{
public:
  Correlation_core_filterbank();
  virtual ~Correlation_core_filterbank();
  virtual void do_task();
  virtual void set_parameters(const Correlation_parameters &parameters,
                              std::vector<Delay_table_akima> &delays,
                              std::vector<std::vector<double> > &uvw,
                              int node_nr, double DM_, 
                              bool no_intra_channel_dedispersion_);
protected:
  void integration_initialise();
  void integration_step(std::vector<Complex_buffer> &integration_buffer, 
                        int index);
  void create_baselines(const Correlation_parameters &parameters);
  void sub_integration();
  void integration_write_subints(std::vector<Complex_buffer> &integration_buffer);
  void create_channel_offsets();
  std::vector<double> offsets;
  int nr_subints_per_integration;
  bool no_intra_channel_dedispersion;
  double DM;
};

#endif /*CORRELATION_CORE_H_*/
