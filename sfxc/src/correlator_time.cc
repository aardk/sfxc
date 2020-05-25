#include <stdio.h>
#include <limits>
#include <stdexcept>
#include "correlator_time.h"

Time::Time() : clock_rate(MAX_SAMPLE_RATE), sample_rate(1), nticks(0){
}

Time::Time(int mjd, double sec, double sample_rate_) : clock_rate(MAX_SAMPLE_RATE){
  sample_rate = (int) round(clock_rate / sample_rate_);
  set_time(mjd, sec);
}

Time::Time(double usec, double sample_rate_) : clock_rate(MAX_SAMPLE_RATE){
  sample_rate = (int) round(clock_rate / sample_rate_);
  nticks = (int64_t)round(usec * (clock_rate / 1000000));
}

Time::Time(const std::string &time) : clock_rate(MAX_SAMPLE_RATE), sample_rate(1){
  struct tm tm;
  const char *end;

  tm.tm_hour = tm.tm_min = tm.tm_sec = 0;
  end = strptime(time.c_str(), "%Yy%jd%Hh%Mm%S", &tm);
  // Skip subsecond field for now
  if (end && *end == '.') {
    end++;
    while (isdigit(*end))
      end++;
  }
  if (end && *end == 's')
    end++;
  if (end == NULL)
    end = strptime(time.c_str(), "%Yy%jd%Hh%Mm", &tm);
  if (end == NULL)
    end = strptime(time.c_str(), "%Yy%jd%Hh", &tm);
  if (end == NULL)
    end = strptime(time.c_str(), "%Yy%jd", &tm);
  if (end == NULL || *end)
    throw std::invalid_argument("Invalid datetime string: " + time);
  int time_mjd = mjd(1, 1, 1900 + tm.tm_year) + tm.tm_yday;
  set_time(time_mjd, 60 * (60 * tm.tm_hour + tm.tm_min) + tm.tm_sec);
}

Time Time::now() {
  struct timeval tv;
  gettimeofday(&tv, NULL);

  struct tm *t = localtime(&tv.tv_sec);
  int m = mjd(t->tm_mday, t->tm_mon + 1, t->tm_year + 1900);
  int sec = 60 * (60 * t->tm_hour + t->tm_min) + t->tm_sec;
  Time tnow(m, sec);
  tnow.inc_time_usec(tv.tv_usec);
  return tnow;
}

void Time::set_sample_rate(double sample_rate_){
  sample_rate = (int) round(clock_rate / sample_rate_);
}

double Time::get_mjd() const{
  int64_t ticks_per_day = (int64_t) (SECONDS_PER_DAY  * clock_rate);
  int ndays = nticks / ticks_per_day;
  double mjd = REFERENCE_MJD + ndays;
  return mjd + (nticks % ticks_per_day) / clock_rate / SECONDS_PER_DAY;
}

double Time::get_time_usec() const{
  int64_t ticks_per_day = (int64_t) (SECONDS_PER_DAY * clock_rate);
  return (nticks % ticks_per_day) / (clock_rate / 1000000);
}

double Time::get_time() const{
  int64_t ticks_per_day = (int64_t) (SECONDS_PER_DAY * clock_rate);
  return (nticks % ticks_per_day) / clock_rate;
}

const Time Time::max_time(){
  Time new_time;
  new_time.nticks = std::numeric_limits<int64_t>::max();
  return new_time;
}

double Time::diff(const Time &other) const {
  return (nticks - other.nticks) / clock_rate;
}

void Time::get_date(int &year, int &day) const{
// Valid from 1-1-1901 to 1-1-2100
  const int mjd1901 = 15385;
  int ndays = get_mjd() - mjd1901;
  int p1 = ndays / (365 * 4 + 1);
  int p2 = (ndays - p1 * (365 * 4 + 1)) / 365;
  year = 4 * p1 + p2 + 1901;
  day = 1 + (ndays - p1 * (365 * 4 + 1)) % 365;
}

void Time::get_time(int &h, int &m, double &s) const{
  double time = (nticks % ((int64_t)clock_rate * 24*60*60)) / clock_rate;
  h = (int) floor(time / (60 * 60));
  time -= h * 60 * 60 ;
  m = (int) floor(time / 60);
  time -= m * 60;
  s = time;
}

std::string Time::date_string(int precision /* = 3 */ ) const{
  int year, day, hour, minute;
  double second;
  char date[40], fmt[40];
  if (precision > 17)
    precision = 17;
  int nlead = (precision > 0) ? precision + 3 : 2;
  snprintf(fmt, 40, "%%04dy%%03dd%%02dh%%02dm%%0%d.%dfs", 
           nlead, precision);

  get_date(year, day);
  get_time(hour, minute, second);
  snprintf(date, 40, fmt,
           year, day, hour, minute, second);
  return std::string(date);
}

std::ostream& operator<<(std::ostream &o, const Time t){
  o << t.date_string();
  return o;
}
