
#include <cstdlib>

#include <libcdsBasics.h>
#include <BitSequence.h>
#include <Mapper.h>
#include <Sequence.h>

using namespace std;
using namespace cds_static;


void testIntersectionWTNoptrs(Array & a,WaveletTreeNoptrs& wt) {
    uint tests = 500;

    while(tests) {
        /* create ranges */
        std::vector< pair<size_t,size_t> > ranges;
        uint N = 2 + rand() % 5;
        for(size_t i=0;i<N;i++) {
            size_t left = rand() % (a.getLength()-10);
            size_t max = a.getLength() - left - 1;
            size_t right = left + 1+ ((rand()%64) % max);
            ranges.push_back( std::make_pair(left,right) );
        }
        /* create a copy */
        std::vector< pair<size_t,size_t> > r2 = ranges;

        /* perform intersection */
        std::vector<uint>* isect = wt.intersect(ranges);

        /* check result */
        for(size_t j=0;j<isect->size();j++) {
            uint res = (*isect)[j];
            uint found = 0;
            for(size_t i=0;i<N;i++) {
                pair<size_t,size_t> cr = r2[i];
                
                for(size_t k=cr.first;k<=cr.second;k++) {
                    if( a[k] == res ) {
                        found++;
                        break;
                    }
                }
            }
            if(found!=N) {
                fprintf(stdout,"INTERSECTION ERROR!\n");
                exit(EXIT_FAILURE);
            }
        }
        delete isect;

        tests--;
    }
}

void testIntersectionWT(Array & a,WaveletTree& wt) {
    uint tests = 500;

    while(tests) {
        /* create ranges */
        std::vector< pair<size_t,size_t> > ranges;
        uint N = 2 + rand() % 5;
        for(size_t i=0;i<N;i++) {
            size_t left = rand() % (a.getLength()-10);
            size_t max = a.getLength() - left - 1;
            size_t right = left + 1+ ((rand()%64) % max);
            ranges.push_back( std::make_pair(left,right) );
        }
        /* create a copy */
        std::vector< pair<size_t,size_t> > r2 = ranges;

        /* perform intersection */
        std::vector<uint>* isect = wt.intersect(ranges);

        /* check result */
        for(size_t j=0;j<isect->size();j++) {
            uint res = (*isect)[j];
            uint found = 0;
            for(size_t i=0;i<N;i++) {
                pair<size_t,size_t> cr = r2[i];
                
                for(size_t k=cr.first;k<=cr.second;k++) {
                    if( a[k] == res ) {
                        found++;
                        break;
                    }
                }
            }
            if(found!=N) {
                fprintf(stdout,"ERROR! could not find %u in enough ranges\n",res);
                exit(EXIT_FAILURE);
            }
        }
        delete isect;

        tests--;
    }
}


int main(int argc, char ** argv) {

  srand( getpid() );

  uint len = 250 + rand() % 50000;
  uint maxv = 20 + rand() % 5;

  Array a(len,maxv);
  for(uint i=0;i<len;i++) {
    a.setField(i,rand()%maxv);
  }
  
  Mapper * mapper = new MapperCont(a, BitSequenceBuilderRG(20));
  Mapper * mapper2 = new MapperNone();
  mapper->use();
  mapper2->use();

  cout << "Test 1 : Wavelet tree without pointers" << endl;
  WaveletTreeNoptrs wt3(a, new BitSequenceBuilderRRR(32), mapper);
  cout << "bs.size() = " << wt3.getSize() << endl;
  testIntersectionWTNoptrs(a,wt3);

  cout << "Test 2 : wavelet tree with pointers" << endl;
  WaveletTree wt1(a,new wt_coder_huff(a,mapper),new BitSequenceBuilderRG(20),mapper);
  cout << "bs.size() = " << wt1.getSize() << endl;
  testIntersectionWT(a,wt1);
  mapper->unuse();
  mapper2->unuse();
  fprintf(stdout,"ALL OK\n");
  return 0;
}

