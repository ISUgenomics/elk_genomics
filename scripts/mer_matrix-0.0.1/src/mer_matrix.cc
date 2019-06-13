#include <config.h>
#include <src/mer_matrix_cmdline.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/stream_iterator.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jflib/multiplexed_io.hpp>

mer_matrix_cmdline args;

using                                                     jflib::o_multiplexer;
using                                                     jellyfish::mer_dna;
typedef std::vector<const char*>                          file_vector;
typedef file_vector::const_iterator                       file_iterator;
typedef jellyfish::stream_iterator<file_iterator>         stream_iterator;
typedef jellyfish::whole_sequence_parser<stream_iterator> sequence_parser;

class compute_matrix : public jellyfish::thread_exec {
  sequence_parser parser_;
  o_multiplexer   multiplexer_;
public:
  compute_matrix(uint16_t nb_threads, const file_vector& paths) :
    parser_(3 * nb_threads, stream_iterator(paths.cbegin(), paths.cend()),
            stream_iterator()),
    multiplexer_(&std::cout, 3 * nb_threads, 4096)
  { }

  virtual void start(int thid) {
    const int nb_mers = 1 << (2 * mer_dna::k());
    jflib::omstream out(multiplexer_);
    uint64_t        counts[nb_mers];
    mer_dna         m;

    while(true) {
      sequence_parser::job j(parser_);
      if(j.is_empty())
        break;
      memset(counts, '\0', sizeof(counts));
      unsigned int filled = 0;
      for(auto base = j->seq.cbegin(); base != j->seq.cend(); ++base) {
        if(m.shift_left(*base) != 'N') {
          if(++filled >= mer_dna::k())
            ++counts[m.word(0)];
        } else {
          filled = 0;
        }
      }
      out << j->header;
      for(int i = 0; i < nb_mers; ++i) {
        out << " " << counts[i];
      }
      out << "\n";
      out << jflib::endr;
    }
  }
};

int main(int argc, char *argv[])
{
  std::ios::sync_with_stdio(false);

  args.parse(argc, argv);
  if(args.mer_arg > 32)
    die << "Mer size limited to 32";
  mer_dna::k(args.mer_arg);

  compute_matrix compute(args.threads_arg, args.sequence_arg);
  compute.exec_join(args.threads_arg);

  return 0;
}
