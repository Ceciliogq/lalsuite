/*  <lalVerbatim file="LALDeterminant3CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALDeterminant.c}}

Module to calculate the determinant of a 3-dimensional matrix $g_{ij}$. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALDeterminant3CP}
\index{\verb&LALDeterminant()& 3D determinant}

\subsubsection*{Description}

This code computes the determinant of a 3-dimensional matrix.

\subsubsection*{Algorithm}
Given a matrix $g_{ij}$ its determinant
is computed using  the formula $g = \epsilon^{ijk} g_{i1} g_{j2} g_{k3},$ 
where $\epsilon$ is the totally anti-symmetric tensor in 3-dimensions. 

\subsubsection*{Uses}
None.

\subsubsection*{Notes}
Don't ever generalise this to higher dimensions since this 
would take many more operations than some of the standard routines.

\vfill{\footnotesize\input{LALDeterminant3CV}}

</lalLaTeX>  */

#include <lal/LALInspiralBank.h>

NRCSID(LALDETERMINANT3C, "$Id$");

/*  <lalVerbatim file="LALDeterminant3CP"> */

void LALDeterminant3(LALStatus *status, 
                     REAL8     *determinant, 
                     REAL8     **matrix) 
{ /* </lalVerbatim> */

   REAL8 epsilon[3][3][3] = {{
                              { 0, 0, 0},
                              { 0, 0, 1},
                              { 0,-1, 0}},
                             {{ 0, 0,-1},
                              { 0, 0, 0},
                              { 1, 0, 0}},
                             {{ 0, 1, 0},
                              {-1, 0, 0},
                              { 0, 0, 0}}};
   INT4 Dim=3,i,j,k;
   INITSTATUS(status, "LALDeterminant3", LALDETERMINANT3C);
   ATTATCHSTATUSPTR(status);

   ASSERT (matrix,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   *determinant = 0;
   for (i=0; i<Dim; i++) {
   for (j=0; j<Dim; j++) {
   for (k=0; k<Dim; k++) {
      *determinant+=epsilon[i][j][k]*matrix[0][i]*matrix[1][j]*matrix[2][k];
   }}}

   ASSERT (*determinant != 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
