#include <Rcpp.h>
#include <set>
using namespace Rcpp;


// [[Rcpp::export]]
int filter_check_dosing_table(DataFrame& dosingTable, NumericVector inputs)
{
  //Pointers to required columns
  double *inputCol = REAL(dosingTable["INPUT"]);
  double *dosingTableTimes = REAL(dosingTable["TIME"]);
  double *dosingTableDuration = REAL(dosingTable["DURATION"]);
  double *dosingTableLagtime = REAL(dosingTable["LAGTIME"]);

  //Iterate for number of dosing inputs times
  for(int i=0; i<inputs.length();i++)
  {
    std::set<double> dosetimes;
    std::vector<double> dosingInputkTime;
    std::vector<double> dosingInputkTestTime;

    // For each k dosing input create lists of Dosing Test Time which is equal to
    // doestime + duration + lagtime
      for(int j=0; j<dosingTable.nrow();j++)
      {
        if(inputCol[j]==inputs[i])
        {
          dosetimes.insert(dosingTableTimes[j]);
          dosingInputkTestTime.push_back(dosingTableTimes[j]+dosingTableDuration[j]+dosingTableLagtime[j]);
          dosingInputkTime.push_back(dosingTableTimes[j]);
        }
      }

      // if dose is present for the current input check if dose administration of a dose happens before start of next dosing event
      // if it happens after start of next dose, show error message and exit
      if(dosetimes.size()>1)
      {
        std::vector<double> inputkTestTimeFiltered;
        for(std::set<double>::iterator it=dosetimes.begin(); it!= --dosetimes.end();it++)
        {
          std::vector<double>::iterator findIterator=dosingInputkTime.begin();
          while(findIterator!=dosingInputkTime.end())
          {
            findIterator = std::find(findIterator,dosingInputkTime.end(),*it);
            if(findIterator==dosingInputkTime.end())
              break;
            inputkTestTimeFiltered.push_back(dosingInputkTestTime[findIterator-dosingInputkTime.begin()]);
            std::advance(findIterator,1);
          }
          std::set<double>::iterator itnext=it;
          std::advance(itnext,1);
          if(*(std::max_element(inputkTestTimeFiltered.begin(),inputkTestTimeFiltered.end()))>*(itnext))
          {
            inputkTestTimeFiltered.clear();
            return 1;
          }
        }
        inputkTestTimeFiltered.clear();
        std::vector<double>().swap(inputkTestTimeFiltered);
      }
      dosingInputkTime.clear();
      std::vector<double>().swap(dosingInputkTime);
      dosingInputkTestTime.clear();
      std::vector<double>().swap(dosingInputkTestTime);
      dosetimes.clear();
      std::set<double>().swap(dosetimes);
  }
  return 0;
}
