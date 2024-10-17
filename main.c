/* main.c
 *
 * Copyright 2024 Stefan Thiago
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define f11 0.9808f
#define f13 -0.008f
#define b11 0.008f
#define b33 0.000012f
#define q1 0.0324f
#define q2 0.0269f
#define q3 0.0211f
#define r1 0.068f

void compute_p(void)
{
  int j, l, k, count = 0;
  float p[4][4] = {
    {1., 1., 1., 1.},
    {1., 1., 1., 1.},
    {1., 1., 1., 1.},
    {1., 1., 1., 1.}};
  float pkb[4][4];
  float kk[4][2];
  float p_old[4][4];
  float aux_m[2][2];
  float aux_p[4][4];
  float aux;
  float diff_p;
  do
  {
    // p_old = p
    for(j=0;j<4;j++)
    {
      for(l=0;l<4;l++)
        p_old[j][l] = p[j][l];
    }

    // pkb = f*p*f' + q
    pkb[0][0] = f11*f11*p[0][0] + f11*f13*(p[0][2]+p[2][0]) + f13*f13*p[2][2];
    pkb[1][0] = f11*f11*p[1][0] + f11*f13*(p[3][0]+p[1][2]) + f13*f13*p[3][2];
    pkb[2][0] = f11*p[2][0] + f13*p[2][2];
    pkb[3][0] = f11*p[3][0] + f13*p[3][2];
    pkb[0][1] = f11*f11*p[0][1] + f11*f13*(p[0][3]+p[2][1]) + f13*f13*p[2][3];
    pkb[1][1] = f11*f11*p[1][1] + f11*f13*(p[3][1]+p[1][3]) + f13*f13*p[3][3];
    pkb[2][1] = f11*p[2][1] + f13*p[2][3];
    pkb[3][1] = f11*p[3][1] + f13*p[3][3];
    pkb[0][2] = f11*p[0][2] + f13*p[2][2];
    pkb[1][2] = f11*p[1][2] + f13*p[3][2];
    pkb[2][2] = p[2][2];
    pkb[3][2] = p[3][2];
    pkb[0][3] = f11*p[0][3] + f13*p[2][3];
    pkb[1][3] = f11*p[1][3] + f13*p[3][3];
    pkb[2][3] = p[2][3];
    pkb[3][3] = p[3][3];

    pkb[0][0] = pkb[0][0] + q1;
    pkb[1][1] = pkb[1][1] + q1;
    pkb[2][2] = pkb[2][2] + q2;
    pkb[3][3] = pkb[3][3] + q3;

    // kk = pkb*h'/(h*pkb*h' + r);
    aux_m[0][0] = pkb[0][0] + r1;
    aux_m[0][1] = pkb[0][1];
    aux_m[1][0] = pkb[1][0];
    aux_m[1][1] = pkb[1][1] + r1;
    aux = 1/(aux_m[0][0]*aux_m[1][1] - aux_m[0][1]*aux_m[1][0]);

    kk[0][0] = aux*(pkb[0][0]*aux_m[1][1] - pkb[0][1]*aux_m[1][0]);
    kk[1][0] = aux*(pkb[1][0]*aux_m[1][1] - pkb[1][1]*aux_m[1][0]);
    kk[2][0] = aux*(pkb[2][0]*aux_m[1][1] - pkb[2][1]*aux_m[1][0]);
    kk[3][0] = aux*(pkb[3][0]*aux_m[1][1] - pkb[3][1]*aux_m[1][0]);
    kk[0][1] = aux*(pkb[0][1]*aux_m[0][0] - pkb[0][0]*aux_m[0][1]);
    kk[1][1] = aux*(pkb[1][1]*aux_m[0][0] - pkb[1][0]*aux_m[0][1]);
    kk[2][1] = aux*(pkb[2][1]*aux_m[0][0] - pkb[2][0]*aux_m[0][1]);
    kk[3][1] = aux*(pkb[3][1]*aux_m[0][0] - pkb[3][0]*aux_m[0][1]);

    // pk = (eye(length(xk)) - kk*h)*pkb
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        aux_p[j][k] = 0;
        for (l = 0; l < 2; l++)
          aux_p[j][k] += kk[j][l] * pkb[l][k];
      }
    }
    diff_p = 0;
    for(j=0;j<4;j++)
    {
      for(l=0;l<4;l++)
      {
        p[j][l] = pkb[j][l] - aux_p[j][l];
        diff_p += fabs(p[j][l] - p_old[j][l]);
      }
    }
    count++;
  }while(diff_p > 1E-6);

  printf("Resultant Matrix is:\n");

  for(j=0;j<4;j++)
  {
    for(l=0;l<4;l++)
      printf("%.3f\t", p[j][l]);
    printf("\n");
  }
  printf("Diff: %f\n",diff_p);
  printf("Count: %d\n", count);
}

int main (int argc, char *argv[])
{
  printf ("Hello, World!\n");
  compute_p();
  return EXIT_SUCCESS;
}
