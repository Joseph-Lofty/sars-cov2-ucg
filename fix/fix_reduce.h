/* ----------------------------------------------------------------------
    The code is for extra fix style developed by Fengyu Li, Yuwei Zhang,
  Fei Xia, and Xin Xu. One can use this fix style to realize anisotropic
  network model. The computation costs can be reduced via a technique of
  grouping pairwise interactions with the same k (force constant) and r0
  (equilibrium distance) with the help of two additional files.
    1. file1 includes the group indices.
      data structure: atom1_index, atom2_index, group_index.
    2. file2 includes specific force field parameters of each group.
      data structure: group_index, k, r0.

                              fxia@chem.ecnu.edu.cn, xxchem@fudan.edu.cn

                                                               Feb. 2023
------------------------------------------------------------------------- */

/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(reduce, FixREDUCE)

#else

#ifndef LMP_FIX_REDUCE_H
#define LMP_FIX_REDUCE_H

#include "fix.h"
#include <string>
using namespace std;

namespace LAMMPS_NS {

class FixREDUCE : public Fix {
public:
  FixREDUCE(class LAMMPS *, int, char **);
  virtual ~FixREDUCE() {
    delete [] list1;
    delete [] list2;
    delete []idx_local;
  }
  int setmask();
  virtual void init();
  virtual void post_force(int);
  void setup(int);
protected:
  struct file1 {
    int atom1_ID;
    int atom2_ID;
    int bond_type;
  };
  struct file2 {
    double k;
    double r0;
  };

  int getNumberOfLine(string);
  void zero_list1(struct file1 * list1);
  void zero_list2(struct file2 * list2);
  void get_list1(struct file1 *, string);
  void get_list2(struct file2 *, string);

  string file1_name;
  string file2_name;
  int natoms;
  int numberOfLine1;
  int numberOfLine2;
  file1 * list1;
  file2 * list2;
  tagint *idx_local;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
