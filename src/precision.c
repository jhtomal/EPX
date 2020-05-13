/*

   precision.c

   Originally written by Marcia Wang, University of Waterloo; general modifications and changing the sorting mechanism using https://phoxis.org/2012/07/12/get-sorted-index-orderting-of-an-array/ by Grace Hsu.

*/



/*

When there are ties in estimated probability of relevance, we find the expected AHR instead of AHR. The algorithm below give us the analytic result.

*/

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
// #define M 7; // for indexx

/* holds the address of the array of which the sorted index
 * order needs to be found
 */
double *base_arr;


/*
void preci();

int main() {
  int n = 3;
  double score[3] = {0.6, 1, 0.9};
  int trueclass[3] = {0, 1, 1};
  double precision = 0;

  preci(&n, score, trueclass, &precision);

  printf("%g\n", precision);
}
*/

double
expection (int i, int j, int n,
           int active[], int total[],
           int group[]);

// void
// indexx(int n, double arr[], int indx[]);

void
sortindex (int n, double arr[], int idx[]);

void
preci(int  *n, double *score,
      int  *trueclass, double *precision )

{

   int nobs= *n, class[nobs], index[nobs],i, j, p,

                  group[nobs], active[nobs], total[nobs];

   double rank[nobs], temp, sum , actnum=0.0;

   //indexx(nobs, score, index);
   sortindex(nobs, score, index);

   /* reorder the response and their rank in an acending order */

   for ( i=0; i< nobs; i++){

      rank[i]=score[index[nobs-i-1]];
      class[i]=trueclass[index[nobs-i-1]];
   }


   p=0;
   group[0]=0;
   active[p]=class[0];
   total[p]=1;

   /*

    "group[i]" will record the group index of the ith obs. The permutation
    will occur within each group (total group number=p, from 0 to p-1). The
    ith obs will be in group[i]th group. so group[1]=0 and group[n-1]=p-1.
    we are going to compute the "group".

   */

    for (i=1; i<nobs; i++){
      if  (fabs(rank[i]-rank[i-1])>1e-7){

        p++;
        active[p]=class[i];
        total[p]=1;

      }

     else {

        active[p]=active[p]+class[i];
        total[p]=total[p]+1;

     }

     group[i]=p;

   }



   for (i=0; i<p+1; i++)

        actnum += (double) (active[i]);

        /*  printf(" %f  \n ", actnum);


            for (i=0; i<nobs; i++)
            printf(" %d \n ", group[i]);
        */
  sum=0.0;

  for (i=0; i<nobs; i++){

     if (active[group[i]] >0){

        temp=0.0;

        for (j=0; j<=i; j++)

           temp += expection(i, j, nobs, active, total, group);


        sum += temp/((double) (i+1));

     }
  }

  /* printf("sum=  %f  sum= \n ", sum, sum/actnum); */

  *precision=sum/actnum;


}


double
expection (int i, int j, int n,

           int active[], int total[],

           int group[])

{

    double expt;

    if (group[i]!=0){

       if (active[group[i]]==1){

         if (group[i]==group[j]){

            if (i==j)

               expt= 1.0/(double) (total[group[i]]);

            else

               expt= 0.0;
         }

         else {

              if( active[group[j]]==0)

                  expt=0.0;

              else

                  expt= ((double) (active[group[i]])/(double) (total[group[i]]))*
                        ((double) (active[group[j]])/(double) (total[group[j]]));
       }

     }

     else {

        if (group[i]==group[j]){

           if (i==j)

              expt= (double) (active[group[i]])/(double) (total[group[i]]);

           else

               expt=  ((double) (active[group[i]])/(double) (total[group[i]]))*
                ((double) (active[group[i]]-1)/(double) (total[group[i]]-1));

        }

        else {

              if( active[group[j]]==0)

                    expt=0.0;

              else

                    expt= ((double) (active[group[i]])/(double) (total[group[i]]))*
                     ((double) (active[group[j]])/(double) (total[group[j]]));

        }

    }

  }

  else {

     if (active[group[i]]==1){

        if (i==j)

           expt= 1.0/(double) (total[group[i]]);

        else

           expt= 0;

     }

     else {

          if (i==j)

            expt= (double) (active[group[i]])/(double) (total[group[i]]);

          else

             expt= ((double) (active[group[i]])/(double) (total[group[i]]))*
               (double) (active[group[i]]-1)/(double) (total[group[i]]-1);

      }

  }

  return(expt);

}

/* CODE BELOW FROM https://phoxis.org/2012/07/12/get-sorted-index-orderting-of-an-array/
*/

/*
Comparison function for sorting compares the value at the original arr indexed by idx. The qsort function calls the comparison function using the values in the passed array to it, which is the pointer to idx.
*/

/* Note how the compare function compares the values of the
 * array to be sorted. The passed value to this function
 * by `qsort' are actually the `idx' array elements.
 */
  /*
  Note that the ugly looking *((int *)a) is nothing but typecasting the void * pointer to integer first, then dereferencing it to get the value pointed to by the pointer.
  */
static int
compar (const void *a, const void *b)
{
  int aa = *((int *) a), bb = *((int *) b);
  if (base_arr[aa] < base_arr[bb])
    return -1;
  if (base_arr[aa] == base_arr[bb])
    return 0;
  if (base_arr[aa] > base_arr[bb])
    return 1;
  else{
    return 0;
  }
}

/*
Here we compare the values in arr, the original array, indexed with the passed integer values in compar function, which are actually the values of the idx array. There is one issue. Where will we get the arr inside the scope of this function? One way is we can make a global or static global pointer pointing to the original array and use it in the compare function.
*/

void
sortindex (int n, double arr[], int idx[])
{
  int i;

  /* initialize initial index permutation of unmodified `arr'
   */
  for (i = 0; i < n; i++)
    {
      idx[i] = i;
    }

  /* Assign the address of out original array to the static global
   * pointer, this will be used by the compare function to index
   * into the original array using `idx' values
   */
  base_arr = arr;

  qsort (idx, n, sizeof (int), compar);


}