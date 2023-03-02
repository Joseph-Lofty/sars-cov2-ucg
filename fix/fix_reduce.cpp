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

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_reduce.h"
#include <cstring>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

#include "comm.h"
#include <fstream>
#include "mpi.h"
#include <cmath>
#include <string>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

int FixREDUCE::getNumberOfLine(string file_name) {
  ifstream ifs;
  ifs.open(file_name, std::ios::in);
  int numberOfLine = 0;
  char c;
  while ((c = ifs.get()) != EOF) {
    if (c == '\n') {
      numberOfLine++;
    }
  }
  ifs.close();
  return numberOfLine;
}

void FixREDUCE::zero_list1(struct file1* list1) {
  for (int i = 0; i < sizeof(list1) / sizeof(list1[0]); i++) {
    list1->atom1_ID = 0;
    list1->atom2_ID = 0;
    list1->bond_type = 0;
  }
}

void FixREDUCE::zero_list2(struct file2* list2) {
  for (int i = 0; i < sizeof(list2) / sizeof(list2[0]); i++) {
    list2->k = 0;
    list2->r0 = 0;
  }
}

void FixREDUCE::get_list1(struct file1* list1, string file1_name) {
  ifstream ifs;
  ifs.open(file1_name, std::ios::in);
  char buff[1024] = { 0 };
  int count = 0;
  while (ifs.getline(buff, sizeof(buff))) {
    sscanf(buff, "%d %d %d", &list1[count].atom1_ID, &list1[count].atom2_ID, &list1[count].bond_type);
    count++;
  }
  ifs.close();
}

void FixREDUCE::get_list2(struct file2* list2, string file2_name) {
  ifstream ifs;
  ifs.open(file2_name, std::ios::in);
  char buff[1024] = { 0 };
  int count = 0;
  while (ifs.getline(buff, sizeof(buff))) {
    sscanf(buff, "%*d %lf %lf", &list2[count].k, &list2[count].r0);
    count++;
  }
  ifs.close();
}

/* ---------------------------------------------------------------------- */

FixREDUCE::FixREDUCE(LAMMPS* lmp, int narg, char** arg) :
  Fix(lmp, narg, arg) {
  if (narg != 5) {
    error->all(FLERR, "\nUsage:\n  fix     fix_name all reduce bond_type_file1_name bond_parm_file2_name\n");
  }
  file1_name = arg[3];
  file2_name = arg[4];
}

/* ---------------------------------------------------------------------- */

int FixREDUCE::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixREDUCE::init() {

  numberOfLine1 = getNumberOfLine(file1_name);
  list1 = new file1[numberOfLine1];
  zero_list1(list1);
  get_list1(list1, file1_name);

  numberOfLine2 = getNumberOfLine(file2_name);
  list2 = new file2[numberOfLine2];
  zero_list2(list2);
  get_list2(list2, file2_name);

  natoms = atom->natoms;
  idx_local = new tagint[natoms];

  if (strstr(update->integrate_style, "respa")) {
    error->all(FLERR, "Fix Reduce does not support respa");
  }
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixREDUCE::post_force(int vflag) {
  // update v and x of atoms in group

  double** x = atom->x;
  double** f = atom->f;
  int* mask = atom->mask;
  int nlocal = atom->nlocal;
  tagint* tag = atom->tag;

  for (int i = 0; i < natoms; i++) {
    idx_local[tag[i] - 1] = i;
  }

  for (int i = 0; i < numberOfLine1; i++) {
    int ii, jj;
    int idx_i, idx_j;

    ii = list1[i].atom1_ID - 1;
    jj = list1[i].atom2_ID - 1;
    idx_i = idx_local[ii];
    idx_j = idx_local[jj];

    if (idx_i < nlocal || idx_j < nlocal) {

      int type = list1[i].bond_type - 1;
      double k = list2[type].k;
      double r0 = list2[type].r0;

      double delx = x[idx_i][0] - x[idx_j][0];
      double dely = x[idx_i][1] - x[idx_j][1];
      double delz = x[idx_i][2] - x[idx_j][2];

      double rsq = delx * delx + dely * dely + delz * delz;
      double r = sqrt(rsq);

      double delr = r - r0;
      double rk = k * delr;

      double fbond;
      if (r > 0.0) fbond = -2.0 * rk / r;
      else fbond = 0.0;

      if (idx_i < nlocal) {
        f[idx_i][0] += delx * fbond;
        f[idx_i][1] += dely * fbond;
        f[idx_i][2] += delz * fbond;
      }

      if (idx_j < nlocal) {
        f[idx_j][0] -= delx * fbond;
        f[idx_j][1] -= dely * fbond;
        f[idx_j][2] -= delz * fbond;
      }

    }
  }
}

void FixREDUCE::setup(int vflag) {
  if (strstr(update->integrate_style, "verlet")) {
    post_force(vflag);
  }
}
