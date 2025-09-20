// radp.hpp

#ifndef RADP_HPP
#define RADP_HPP


#include "orpp/overallriskdp.hpp"
namespace orpp 
{

using critcvar = CVaR<ldistribution<double>,true>;
class critmcv: public MeanCVaR<ldistribution<double>,true>
{
public:
    critmcv(double alpha) : MeanCVaR(0.95, alpha) {}
};

class critmcv75: public MeanCVaR<ldistribution<double>,true>
{
public:
    critmcv75(double alpha) : MeanCVaR(0.75, alpha) {}
};



}
#endif

