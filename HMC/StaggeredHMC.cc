/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/GenericHmcRunner.h

Copyright (C) 2024

Author: Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory

*************************************************************************************/
/*
    @file StaggeredHMC.cc
    @brief Runs generic staggered HMC
    @author Curtis Taylor Peterson
*/
#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;
  
  Grid_init(&argc, &argv);

  //////////////////////
  // Generic typedefs //
  //////////////////////

  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;
  typedef StaggeredImplR FermionImplPolicy;
  typedef NaiveStaggeredFermionD FermionAction;
  typedef typename FermionAction::FermionField FermionField;
  typename FermionAction::ImplParams params;
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;

  ///////////////////////////////////////
  // Construct generic HMC information //
  ///////////////////////////////////////

  Real trajectory_length = 1.0;
  double cg_tol          = 1e-8;
  double max_cg_iters    = 2000;
  int level2_multiplier  = 4;
  int md_steps           = 20;

  IntegratorParameters MD;
  HMCparameters HMCparams;
  CheckpointerParameters CPparams;
  RNGModuleParameters RNGpar;
  HMCWrapper TheHMC;
  ConjugateGradient<FermionField> CG(cg_tol, max_cg_iters);
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(level2_multiplier);

  TheHMC.Parameters.MD.MDsteps = md_steps;
  TheHMC.Parameters.MD.trajL   = trajectory_length;

  TheHMC.Resources.AddFourDimGrid("gauge");

  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 5;
  CPparams.format = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  TheHMC.Resources.AddObservable<PlaqObs>();

  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  ///////////////////////
  // Construct action //
  //////////////////////

  RealD beta    = 7.0;
  RealD mass    = 0.1;
  RealD c1      = 0.0;
  RealD u0      = 1.0; 
  bool smeared  = false;

  SymanzikGaugeActionR GaugeAction(beta);
  LatticeGaugeField U(GridPtr);
  FermionAction Ds(U, *GridPtr, *GridRBPtr, mass, c1, u0, params); 
  TwoFlavourEvenOddPseudoFermionAction<FermionImplPolicy> Nf4(Ds, CG, CG);

  Nf4.is_smeared = smeared; 

  Level1.push_back(&Nf4);  
  Level2.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  /////////////
  // Run HMC //
  /////////////

  //TheHMC.ReadCommandLine(argc, argv);
  TheHMC.Run();

  Grid_finalize();
}
