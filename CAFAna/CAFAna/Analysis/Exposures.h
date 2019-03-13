#pragma once

namespace ana

{

  //Beam power POT / year
  //*======================================================================================================*/
  //*                                     Beam Power Configurations                                        */
  //*                            (Use this chart to detemine @power in units of pot/yr)                    */
  //*  MI prot/pulse    Energy (GeV)    Cycle time    Beam power (MW)    Uptime&efficiency    pot/year     */
  //*     7.50E+13          120            1.2           1.20E+00              0.56           1.10E+21     */
  //*     7.50E+13           80            0.9           1.07E+00              0.56           1.47E+21     */
  //*     7.50E+13           60            0.7           1.03E+00              0.56           1.89E+21     */
  //*                                                                                                      */
  //*======================================================================================================*/

  const double POT120 = 1.1e21;
  const double POT80 = 1.47e21;
  const double POT60 = 1.89e21;
}
