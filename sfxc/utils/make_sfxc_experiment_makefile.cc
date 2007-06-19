#include <iostream>
#include <fstream>
#include <vex_parser/vexplus.h>
#include <json/json.h>
#include <basic.h>
#include <assert.h>

std::vector< std::vector<std::string> > channels;

void print_general_section(std::ostream &out, 
			   char this_program[],
			   char vex_file[], 
                           char exp_ctrl_file[]);
void print_ccf_section(std::ostream &out, 
		       VexPlus const &vex, 
                       Json::Value const &sfxc_ctrl);
void print_del_section(std::ostream &out, 
		       VexPlus const &vex, 
                       Json::Value const &sfxc_ctrl);
void print_correlate_section(std::ostream &out,
                             VexPlus const &vex, 
                             Json::Value const &sfxc_ctrl);

void print_html_section(std::ostream &out,
                        VexPlus const &vex, 
                        Json::Value const &sfxc_ctrl);

int main(int argc, char *argv[]) {
  if ((argc != 3) && (argc != 4)) {
    std::cout << argv[0] << " " << "<vex-file> <sfxc-exp-ctrlfile> [<outdir>]" 
	      << std::endl;
    exit(1);
  }

  char *this_program  = argv[0];
  char *vex_file      = argv[1];
  char *ctrl_file     = argv[2];
  char *outdir        = ".";
  if (argc == 4) outdir = argv[3];
  
  std::string outfile = outdir; outfile.append("/Makefile");
  std::ofstream out(outfile.c_str());

  VexPlus vex;
  Json::Value sfxc_ctrl;

  // Vex file
  vex = VexPlus(vex_file);
  vex.parseVex();

  { // sfxc_control
    Json::Reader reader;
    std::ifstream in(ctrl_file);
    bool ok = reader.parse(in, sfxc_ctrl);
    if ( !ok ) {
      // report to the user the failure and their locations in the document.
      std::cout  << "Failed to parse sfxc-control file\n"
                 << reader.getFormatedErrorMessages();
      return 1;
    }
  }

  channels = get_channels(vex, sfxc_ctrl);

  print_general_section(out, this_program, vex_file, ctrl_file);
  print_ccf_section(out, vex, sfxc_ctrl);
  print_del_section(out, vex, sfxc_ctrl);
  print_correlate_section(out, vex, sfxc_ctrl);
  print_html_section(out, vex, sfxc_ctrl);
}

void print_general_section(std::ostream &out, 
			   char this_program[],
			   char vex_file[], 
                           char exp_ctrl_file[]) {
  out << "# SFXC - makefile for generating the ccf-files and delay tables"
      << std::endl;
  out << "# Makefile generated by " << this_program << std::endl;
  out << std::endl;
  out << "VEXFILE =   " << vex_file << std::endl;
  out << "CTRLFILE =  " << exp_ctrl_file << std::endl;
  out << std::endl;
  out << "MAKE_MAKEFILE =  " << this_program << std::endl;
  out << ".PRECIOUS: Makefile" << std::endl;
  out << std::endl;
  out << "all: Makefile ccf del" << std::endl;
  out << std::endl;
  out << "Makefile: $(VEXFILE) $(CTRLFILE)" << std::endl;
  out << "\t$(MAKE_MAKEFILE) $(VEXFILE) $(CTRLFILE) ." <<std::endl;
  out << std::endl;
}

void print_ccf_section(std::ostream &out,
                       VexPlus const &vex, 
                       Json::Value const &sfxc_ctrl) {
  out << "CCF = ";
  for (size_t channel = 0; channel < channels.size(); channel++) {
    if (channels[channel].size() == 1) {
      out << " \\\n  "
          << generate_ccf_filename(sfxc_ctrl["ccfdir"].asString(),
                                   vex.ExperName(),
                                   sfxc_ctrl["scan"].asString(),
                                   channels[channel][0], "");
    } else {
      out << " \\\n  "
          << generate_ccf_filename(sfxc_ctrl["ccfdir"].asString(),
                                   vex.ExperName(),
                                   sfxc_ctrl["scan"].asString(),
                                   channels[channel][0],
                                   channels[channel][1]);
    }
  }
  out << std::endl << std::endl;
  out << "ccf: $(CCF)" << std::endl;

  for (size_t channel = 0; channel < channels.size(); channel++) {
    if (channels[channel].size() == 1) {
      out << generate_ccf_filename(sfxc_ctrl["ccfdir"].asString(),
                                   vex.ExperName(),
                                   sfxc_ctrl["scan"].asString(),
                                   channels[channel][0], "");
    } else {
      out << generate_ccf_filename(sfxc_ctrl["ccfdir"].asString(),
                                   vex.ExperName(),
                                   sfxc_ctrl["scan"].asString(),
                                   channels[channel][0],
                                   channels[channel][1]);
    }
    out << ": $(VEXFILE) $(CTRLFILE)" << std::endl
        << "\tsfxc_vex2ccf $(VEXFILE) $(CTRLFILE)"
        << std::endl
        << std::endl;
  }
}

void print_del_section(std::ostream &out,
                       VexPlus const &vex, 
                       Json::Value const &sfxc_ctrl) {

  out << "DELAYTABLES =";
  for (size_t channel = 0; channel < channels.size(); channel++) {
//     for (int station=0; station<vex.N_Stations(); station++) {
      out << " \\\n  "
          <<generate_del_filename(sfxc_ctrl["deldir"].asString(),
                                  vex.ExperName(),
                                  sfxc_ctrl["scan"].asString(),
                                  channels[channel][0],
                                  vex.Station(vex.N_Stations()-1)).c_str();
//     }
  }
  out << std::endl << std::endl;
  out << "del: $(DELAYTABLES)" << std::endl;

  for (size_t channel = 0; channel < channels.size(); channel++) {
    assert(sfxc_ctrl["stations"].isArray());
    out <<generate_dcf_filename(sfxc_ctrl["deldir"].asString(),
                                vex.ExperName(),
                                sfxc_ctrl["scan"].asString(),
                                channels[channel][0]).c_str();
    out << ": $(VEXFILE)" << std::endl
        << "\tsfxc_vex2dcf $(VEXFILE) $(CTRLFILE)"
        << std::endl;

    for (int station=0; station<vex.N_Stations(); station++) {
      out <<generate_del_filename(sfxc_ctrl["deldir"].asString(),
                                  vex.ExperName(),
                                  sfxc_ctrl["scan"].asString(),
                                  channels[channel][0],
                                  vex.Station(station)).c_str()
          << ": " 
          << generate_dcf_filename(sfxc_ctrl["deldir"].asString(),
                                   vex.ExperName(),
                                   sfxc_ctrl["scan"].asString(),
                                   channels[channel][0]).c_str()
          << std::endl
          << "\tdelmo "
          << generate_dcf_filename(sfxc_ctrl["deldir"].asString(),
                                   vex.ExperName(),
                                   sfxc_ctrl["scan"].asString(),
                                   channels[channel][0]).c_str()
          << " .250"
          << std::endl;
    }
  }
  out << std::endl << std::endl;
}

void print_correlate_section(std::ostream &out,
                             VexPlus const &vex, 
                             Json::Value const &sfxc_ctrl) {
  out << "CORRELATION_OUTPUT_FILES = ";
  for (size_t channel = 0; channel < channels.size(); channel++) {
    if (channels[channel].size() == 2) {
      out << " \\\n  "
          << generate_cor_filename(sfxc_ctrl["outdir"].asString(),
                                   vex.ExperName(),
                                   sfxc_ctrl["scan"].asString(),
                                   channels[channel][0], 
                                   channels[channel][1]).c_str();
    } else {
      // Next channel is not the other polarisation
      out << " \\\n  "
          << generate_cor_filename(sfxc_ctrl["outdir"].asString(),
                                   vex.ExperName(),
                                   sfxc_ctrl["scan"].asString(),
                                   channels[channel][0], 
                                   "").c_str();
    }
  }
  out << std::endl
      << "correlate: all"  << std::endl
      << "\t$(MAKE) do_correlate" << std::endl
      << "do_correlate: $(CORRELATION_OUTPUT_FILES)" << std::endl
      << std::endl;
  
  for (size_t channel = 0; channel < channels.size(); channel++) {
    std::string cor_file, ccf_file;
    if (channels[channel].size() == 2) {
      cor_file = generate_cor_filename(sfxc_ctrl["outdir"].asString(),
                                       vex.ExperName(),
                                       sfxc_ctrl["scan"].asString(),
                                       channels[channel][0], 
                                       channels[channel][1]).c_str();
      ccf_file = generate_ccf_filename(sfxc_ctrl["ccfdir"].asString(),
                                       vex.ExperName(),
                                       sfxc_ctrl["scan"].asString(),
                                       channels[channel][0], 
                                       channels[channel][1]).c_str();
    } else {
      // Next channel is not the other polarisation
      cor_file = generate_cor_filename(sfxc_ctrl["outdir"].asString(),
                                       vex.ExperName(),
                                       sfxc_ctrl["scan"].asString(),
                                       channels[channel][0], 
                                       "");
      ccf_file = generate_ccf_filename(sfxc_ctrl["ccfdir"].asString(),
                                       vex.ExperName(),
                                       sfxc_ctrl["scan"].asString(),
                                       channels[channel][0], 
                                       "");
    }
    out << cor_file << ": " << ccf_file;
    for (int station=0; station<vex.N_Stations(); station++) {
      out << " " << generate_del_filename(sfxc_ctrl["deldir"].asString(),
                                          vex.ExperName(),
                                          sfxc_ctrl["scan"].asString(),
                                          channels[channel][0],
                                          vex.Station(station));
    }
    out << std::endl
        << "\tsfxc_SC " << ccf_file << std::endl;;
  }
}


void print_html_section(std::ostream &out,
                        VexPlus const &vex, 
                        Json::Value const &sfxc_ctrl) {
  out << "html: $(CORRELATION_OUTPUT_FILES)" << std::endl
      << "\tcd " << sfxc_ctrl["htmldir"].asString() 
      << "&& produce_html_plotpage $(CCF)"
      << std::endl << std::endl;
}
