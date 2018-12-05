#include <limits.h>
#include <R.h>
#define z(i,j)   zz[i + nrow * j]

static void edge(int *vertex, int p0, int q0);
static void findclassC(double *zz, int nrow, int ncol, int *vertex);
static int findindexC(int *A, int nA, int a, int *index);
static int neighborC(int *patch, int npatch, int nrow, int ncol, int interior, int *neighbor);
static int allC(double *zz, int *neighbor, int nneighboor, double zpatch0);

static int setdiffC(int *A, int nA, int *B, int nB);
static int uniqueC(int *A, int nA);
static int intersectC(int *A, int nA, int *B, int nB);
static int outerC(int *A, int nA, int *B, int nB, int *OUTER);

void extrema2dC(double *zz, int *nnrow, int *nncol,
              int *maxindex, int *nnmax, int *ttotalmax, int *minindex, int *nnmin, int *ttotalmin) /*, int *vertex) */
{
    int nrow = *nnrow, ncol = *nncol, nvert = nrow * ncol;
    int i, j, nbound=2*(nrow+ncol-2), ninterior, npatch, nneighbor, nmax, nmin, totalmax, totalmin;
    
    /* int nb[8]    = {1, ncol-2,      1, nrow-2,           1,     ncol-2,              1,        nrow-2};  */
    int cumnb[8] = {0,      1, ncol-1,   ncol, ncol+nrow-2, ncol+nrow-1, 2*ncol+nrow-3, 2*ncol+nrow-2};    
 
    int *vertexbound = (int *) R_alloc(nbound, sizeof(int));
            
    vertexbound[cumnb[0]] = 0;   
    for (i=1; i<=(ncol-2); i++) vertexbound[cumnb[0]+i] = nrow * i; 
    vertexbound[cumnb[2]] = nrow * (ncol-1);   
    for (i=1; i<=(nrow-2); i++) vertexbound[cumnb[2]+i] = i + nrow * (ncol - 1); 
    vertexbound[cumnb[4]] = nrow - 1 + nrow * (ncol - 1);   
    for (i=1; i<=(ncol-2); i++) vertexbound[cumnb[4]+i] = nrow - 1 + nrow * i; 
    vertexbound[cumnb[6]] = nrow - 1;   
    for (i=1; i<=(nrow-2); i++) vertexbound[cumnb[6]+i] = i; 
                  
    /**/int *vertex = (int *) R_alloc(nvert, sizeof(int)); /**/
    int *uvertex = (int *) R_alloc(nvert, sizeof(int));
    
    /*  if(vertex != NULL) Rprintf("%c memory for vertex \n", 1);   */
    /*  if(uvertex != NULL) Rprintf("%c memory for unique \n", 1);  */
    /*  Rprintf("%d\n", i); */
            
    /* Finding equivalence class induced by neighborhood relation */ 
    /* Neighborhood is defined if horizontal, vertical, diagonal adjacent pixels are the same. */
    /* Partition the image by Rem's algorithm */
     
    findclassC(zz, nrow, ncol, vertex);

    for (i=0; i<nvert; i++)
        uvertex[i] = vertex[i];
        
    /* Finding neighbors of patches including boundary vertex. */

    for (i=0; i<nbound; i++) 
        vertexbound[i] = vertex[vertexbound[i]];

    nbound = uniqueC(vertexbound, nbound);
    
    nmax = 0;
    nmin = 0;
    totalmax = 0;
    totalmin = 0;
   
    int *patch = (int *) R_alloc(nvert, sizeof(int));
    int *neighbor = (int *) R_alloc(nvert, sizeof(int));
       
    for (i=0; i<nbound; i++) {
        npatch = findindexC(vertex, nvert, vertexbound[i], patch);       
        nneighbor = neighborC(patch, npatch, nrow, ncol, 0, neighbor);
        if (nneighbor > 0) {
            if(allC(zz, neighbor, nneighbor, zz[patch[0]]) == nneighbor) {
                for (j=0; j<npatch; j++) 
                    maxindex[totalmax + j] = patch[j];
                    
                totalmax = totalmax + npatch;
                nmax = nmax + 1;
            }
            else if(allC(zz, neighbor, nneighbor, zz[patch[0]]) == -nneighbor) {
                for (j=0; j<npatch; j++) 
                    minindex[totalmin + j] = patch[j];
                                  
                totalmin = totalmin + npatch;
                nmin = nmin + 1;                
            }
        }  
    }


    ninterior = uniqueC(uvertex, nvert);
    ninterior = setdiffC(uvertex, ninterior, vertexbound, nbound);   
    if (ninterior != 0) {

        for (i=0; i<ninterior; i++) {
            npatch = findindexC(vertex, nvert, uvertex[i], patch);
            nneighbor = neighborC(patch, npatch, nrow, ncol, 1, neighbor);
            if (nneighbor > 0) {
                if(allC(zz, neighbor, nneighbor, zz[patch[0]]) == nneighbor) {
                    for (j=0; j<npatch; j++) 
                        maxindex[totalmax + j] = patch[j];
                        
                    totalmax = totalmax + npatch;
                    nmax = nmax + 1;
                }
                else if(allC(zz, neighbor, nneighbor, zz[patch[0]]) == -nneighbor) {
                    for (j=0; j<npatch; j++) 
                        minindex[totalmin + j] = patch[j];
                        
                    totalmin = totalmin + npatch;
                    nmin = nmin + 1;                
                }
            }
        }   
    }
   
    *nnmax = nmax;
    *nnmin = nmin;
    *ttotalmax = totalmax;
    *ttotalmin = totalmin;
}

void findclassCC(double *zz, int *nrow, int *ncol, int *vertex)
{
    findclassC(zz, *nrow, *ncol, vertex);
}
static void findclassC(double *zz, int nrow, int ncol, int *vertex)
{
    /* zz is the image, which are nrow x ncol matrix. */
    
    /* Vertices with same number vertex[j] are in the same patch.       */
    /* Patches are found using Rem's algorithm for finding equivalence  */
    /* classes (in A Dicipline of Programming).                         */
    
    int nvert = nrow * ncol;  /* number of vertex  */
    int i, j;

    for (j = 0; j < nvert; j++)
        vertex[j] = j;

    for (i = 1; i < nrow; i++)
        if (z(i,0) == z((i-1),0))
            edge(vertex, i, (i-1));

    for (j = 1; j < ncol; j++)
        if (z(0,j) == z(0,(j-1)))
            edge(vertex, nrow * j, nrow * (j-1));
            
    for (i = 1; i < nrow; i++)
        for (j = 1; j < ncol; j++) {
            if (z(i,j) == z((i-1),j))
                edge(vertex, i + nrow * j, i-1 + nrow * j);
            if (z(i,j) == z(i,(j-1)))
                edge(vertex, i + nrow * j, i + nrow * (j-1));
            if (z(i,j) == z((i-1),(j-1)))
                edge(vertex, i + nrow * j, i-1 + nrow * (j-1));
                
            /* ADDING */       
            if (z((i-1),j) == z(i,(j-1)))                             
                edge(vertex, i-1 + nrow * j, i + nrow * (j-1));  
            /* ADDING */
        }

    for (j = 0; j < nvert; j++)
        vertex[j] = vertex[vertex[j]];
}
        
static void edge(int *vertex, int p0, int q0)
{
    int p1, q1;

    p1 = vertex[p0];
    q1 = vertex[q0];

    while (p1 != q1)
        if (q1 < p1) {
            vertex[p0] = q1;
            p0 = p1;
            p1 = vertex[p1];
        } else {
            vertex[q0] = p1;
            q0 = q1;
            q1 = vertex[q1];
        }
}

void findindexCC(int *A, int *nA, int *a, int *index, int *nindex)
{
    *nindex = findindexC(A, *nA, *a, index);
}
static int findindexC(int *A, int nA, int a, int *index)
{

    /* Finding index i such that A[i] = a, which is a given integer.    */
    /* A is the set of integers which is                                */
    /*     (1) nonempty                                                 */
    /*     (2) unsorted or sorted                                       */
    /*     (3) replcates for some elements.                             */
    /* nA is number less than the length of input A.                    */
    /*     One may calculate the result not for all elements of input A.*/
                                                                 
    /* index is the set {i : A[i] = a}.                                 */ 
    /* nindex is the number of the set {i : A[i] = a}.                  */
    
    int i, nindex;

    nindex = 0;
    for (i=0; i<nA; i++) 
        if (A[i] == a) 
        {
            index[nindex] = i;
            nindex = nindex + 1;
        }
        
    return(nindex);  
}
    
void neighborCC(int *patch, int *npatch, int *nrow, int *ncol, int *interior, int *neighbor, int *nneighbor)
{
    *nneighbor=neighborC(patch, *npatch, *nrow, *ncol, *interior, neighbor);
}
static int neighborC(int *patch, int npatch, int nrow, int ncol, int interior, int *neighbor) 
{
    /* Finding neighbor of patch including bounday vertex */
    int i, j, nntotal=0; 
    
    int tmpnpatch = npatch;
    int *tmppatch = (int *) R_alloc(tmpnpatch, sizeof(int));
    
    for (i=0; i<tmpnpatch; i++)
        tmppatch[i] = patch[i];
    
    if (!interior) {
        int nb[8]    = {1, ncol-2,      1, nrow-2,           1,     ncol-2,              1,        nrow-2};  
        int cumnb[8] = {0,      1, ncol-1,   ncol, ncol+nrow-2, ncol+nrow-1, 2*ncol+nrow-3, 2*ncol+nrow-2};    
 
        int *boundary = (int *) R_alloc(2*(ncol+nrow-2), sizeof(int));
                
        boundary[cumnb[0]] = 0;   
        for (i=1; i<=(ncol-2); i++) boundary[cumnb[0]+i] = nrow * i; 
        boundary[cumnb[2]] = nrow * (ncol-1);   
        for (i=1; i<=(nrow-2); i++) boundary[cumnb[2]+i] = i + nrow * (ncol - 1); 
        boundary[cumnb[4]] = nrow - 1 + nrow * (ncol - 1);   
        for (i=1; i<=(ncol-2); i++) boundary[cumnb[4]+i] = nrow - 1 + nrow * i; 
        boundary[cumnb[6]] = nrow - 1;   
        for (i=1; i<=(nrow-2); i++) boundary[cumnb[6]+i] = i;  

        int nn[9]    = {3, 5, 3,  5,  3,  5,  3,  5,  8};        
        int cumnn[9] = {0, 3, 8, 11, 16, 19, 24, 27, 32};
        
        int *neighborindex = (int *) R_alloc(40, sizeof(int));     
                              
                              neighborindex[0] = nrow  ; 
        neighborindex[1] = 1; neighborindex[2] = nrow+1;      
                       
        neighborindex[3] = -nrow  ;                       neighborindex[4] = nrow  ; 
        neighborindex[5] = -nrow+1; neighborindex[6] = 1; neighborindex[7] = nrow+1;
        
        neighborindex[8] = -nrow  ; 
        neighborindex[9] = -nrow+1; neighborindex[10] =  1;
        
        neighborindex[11] = -nrow-1; neighborindex[12] = -1; 
        neighborindex[13] = -nrow  ; 
        neighborindex[14] = -nrow+1; neighborindex[15] = 1;
        
        neighborindex[16] = -nrow-1; neighborindex[17] = -1; 
        neighborindex[18] = -nrow  ;
        
        neighborindex[19] = -nrow-1; neighborindex[20] = -1; neighborindex[21] = nrow-1; 
        neighborindex[22] = -nrow  ;                         neighborindex[23] = nrow  ;
        
        neighborindex[24] = -1; neighborindex[25] = nrow-1; 
                                neighborindex[26] = nrow  ;
                                  
        neighborindex[27] = -1; neighborindex[28] = nrow-1; 
                                neighborindex[29] = nrow  ; 
        neighborindex[30] =  1; neighborindex[31] = nrow+1;
        
        neighborindex[32] = -nrow-1; neighborindex[33] = -1; neighborindex[34] = nrow-1;
        neighborindex[35] = -nrow  ;                         neighborindex[36] = nrow  ;
        neighborindex[37] = -nrow+1; neighborindex[38] =  1; neighborindex[39] = nrow+1;  
        
        for (i = 0; i < 8; i++) { 
            if (tmpnpatch != 0) {            
                nb[i] = intersectC(&boundary[cumnb[i]], nb[i], tmppatch, tmpnpatch);
                if (nb[i] != 0) { 
                    tmpnpatch = setdiffC(tmppatch, tmpnpatch, &boundary[cumnb[i]], nb[i]); 
                    
                    int *tmpneighbor = (int *) calloc(nb[i]*nn[i], sizeof(int));
                    nn[i] = outerC(&boundary[cumnb[i]], nb[i], &neighborindex[cumnn[i]], nn[i], tmpneighbor); 

                    for (j=0; j<nn[i]; j++) 
                        neighbor[nntotal + j] = tmpneighbor[j];
                    free(tmpneighbor);
 
                    nntotal = nntotal + nn[i];                     
                    nntotal = setdiffC(neighbor, nntotal, patch, npatch);
                 
                    if (nntotal != 0) 
                        nntotal = uniqueC(neighbor, nntotal); 
                }
            }
        }
                        
        if (tmpnpatch != 0) {
            int *tmpneighbor = (int *) calloc(tmpnpatch*nn[8], sizeof(int));
            nn[8] = outerC(tmppatch, tmpnpatch, &neighborindex[cumnn[8]], nn[8], tmpneighbor); 
                  
            for (j=0; j<nn[8]; j++)
                neighbor[nntotal + j] = tmpneighbor[j];           
            free(tmpneighbor);
       
            nntotal = nntotal + nn[8];      
            nntotal = setdiffC(neighbor, nntotal, patch, npatch);

            if (nntotal != 0)
                nntotal = uniqueC(neighbor, nntotal);            
        }  
 
      
    } else {      
        int nn = 8;
        int *neighborindex = (int *) R_alloc(nn, sizeof(int));
        
        neighborindex[0] = -nrow-1; neighborindex[1] = -1; neighborindex[2] = nrow-1;
        neighborindex[3] = -nrow;                          neighborindex[4] = nrow;
        neighborindex[5] = -nrow+1; neighborindex[6] = 1;  neighborindex[7] = nrow+1;         
        
        int *tmpneighbor = (int *) calloc(tmpnpatch*nn, sizeof(int));
        nn = outerC(tmppatch, tmpnpatch, neighborindex, nn, tmpneighbor); 
        
        for (j=0; j<nn; j++)
            neighbor[nntotal + j] = tmpneighbor[j];           
        free(tmpneighbor);        
          
        nntotal = uniqueC(neighbor, nn);     
        nntotal = setdiffC(neighbor, nntotal, patch, npatch);
    }    

    return(nntotal); 
}                                                        
    
static int allC(double *zz, int *neighbor, int nneighbor, double zpatch0)
{
    int i, maximum = 0;
    
    for (i=0; i<nneighbor; i++) {
        if (zz[neighbor[i]] < zpatch0) {
            if (i != maximum) break;
            maximum = maximum + 1;
        }
        else if (zz[neighbor[i]] > zpatch0) {
            if (i != -maximum) break;
            maximum = maximum - 1;
        } 
    }     
    return(maximum);
}


void setdiffCC(int *A, int *nA, int *B, int *nB)
{
    *nA = setdiffC(A, *nA, B, *nB);
}
static int setdiffC(int *A, int nA, int *B, int nB)
{

    /* Finding set A - B.                                                       */
    /* A and B are the set of integers which is                                 */
    /*     (1) nonempty                                                         */
    /*     (2) sorted or unsorted                                               */
    /*     (3) no replcates.                                                    */
    /* nA and nB are numbers less than the length of input A and B.             */
    /*     One may calculate the result not for all elements of input A and B.  */
                                                                 
    /* A[0:(nA-1)] is the result set A - B, which overwrites input A.                           */ 
    /*     (1) nonempty or empty                                                                */
    /*     (2) sorted                                                                           */
    /*     (3) no replcates for some elements.                                                  */    
    /*     So, if you want to keep the input A, do not use this function.                       */
    /* nA is the number of the result set A - B, which overwrites the length of input A.        */

    /* Window, INT_MIN is -2^31 = -2147483648.                                                */
        
    int lower = 0, upper, j = 0, i, midpoint;
 
    R_isort(A, nA); 
    R_isort(B, nB); 
     
    nA = nA-1;
    nB = nB-1;
        
    if (nA <= nB) {
        for (i=0; i<=nA; i++) 
        { 
            if (A[i] >= B[lower] && A[i] <= B[nB]) {
                upper = nB;  
                if (upper == lower) 
                {
                    if (A[i] == B[lower]) 
                    {
                        A[i] = INT_MIN;
                        break;
                    }            
                } else if (upper == lower + 1) 
                {
                    if (A[i] == B[lower]) 
                    {
                        A[i] = INT_MIN;
                        lower = lower + 1;
                    } else if (A[i] == B[upper]) 
                    {
                        A[i] = INT_MIN;
                        break;
                    }           
                }
                while ((upper - lower) > 1)
                {
                    midpoint = (lower + upper) >> 1;
                    if (A[i] == B[lower]) 
                    {
                        A[i] = INT_MIN;
                        lower = lower + 1;
                        break;
                    } else if (A[i] == B[upper]) 
                    {
                        A[i] = INT_MIN;
                        lower = upper; 
                        break;
                    } else if (A[i] == B[midpoint]) 
                    {
                        A[i] = INT_MIN;
                        lower = midpoint + 1;
                        break;
                    } else if (A[i] > B[midpoint])
                        lower = midpoint; 
                    else 
                        upper = midpoint;
                }
            }
        }
    } else 
    {
        for (i=0; i<=nB; i++) 
        { 
            if (B[i] >= A[lower] && B[i] <= A[nA]) 
            {
                upper = nA;     
                if (upper == lower) 
                {
                    if (B[i] == A[lower]) 
                    {
                        A[lower] = INT_MIN;
                        break;
                    }            
                } else if (upper == lower + 1) 
                {
                    if (B[i] == A[lower]) 
                    {
                        A[lower] = INT_MIN;
                        lower = lower + 1;
                    } else if (B[i] == A[upper]) 
                    {
                        A[upper] = INT_MIN;
                        break;
                    }           
                }
                while ((upper - lower) > 1)
                {
                    midpoint = (lower + upper) >> 1;
                    if (B[i] == A[lower]) 
                    {
                        A[lower] = INT_MIN;
                        lower = lower + 1;
                        break;
                    } else if (B[i] == A[upper]) 
                    {
                        A[upper] = INT_MIN;
                        lower = upper;
                        break;
                    } else if (B[i] == A[midpoint]) 
                    {
                        A[midpoint] = INT_MIN;
                        lower = midpoint + 1;
                        break;
                    } else if (B[i] > A[midpoint])
                        lower = midpoint; 
                    else 
                        upper = midpoint; 
                }
            }
        }
    }
    
    for (i=0; i<=nA; i++)
        if (A[i] != INT_MIN) 
        {
            A[j] = A[i];
            j = j + 1;
        }   
    nA = j;
    return(nA);    
}    

void uniqueCC(int *A, int *nA)
{
    *nA = uniqueC(A, *nA);
}
static int uniqueC(int *A, int nA)
{
    /* Finding unique elements of A.                                        */
    /* A is the set of integers which is                                    */
    /*     (1) nonempty                                                     */
    /*     (2) unsorted or sorted                                           */
    /*     (3) replcates for some elements.                                 */
    /* nA is number less than the length of input A.                        */
    /*     One may calculate the result not for all elements of input A.    */
                                                                 
    /* A[0:(nA-1)] is the sorted result, which overwrites input A.                              */ 
    /*     (1) nonempty                                                                         */
    /*     (2) sorted                                                                           */
    /*     (3) no replcates for some elements.                                                  */    
    /*     So, if you want to keep the input A, do not use this function.                       */
    /* nA is the number of sorted result A[0:(nA-1)], which overwrites the length of input A.   */
    
    int i, j = 0;

    R_isort(A, nA);  
    if (nA != 1) {
        for(i=0; i<(nA-1); i++) 
            if(A[i] != A[i+1]) {
                A[j] = A[i];
                j = j + 1;
            }
        A[j] = A[nA-1];
        nA = j + 1;            
    }
    return(nA);
}

void intersectCC(int *A, int *nA, int *B, int *nB)
{
    *nA = intersectC(A, *nA, B, *nB);
}
static int intersectC(int *A, int nA, int *B, int nB)
{

    /* Finding set A intersect B.                                               */
    /* A and B are the set of integers which is                                 */
    /*     (1) nonempty                                                         */
    /*     (2) sorted or unsorted                                               */
    /*     (3) no replcates.                                                    */
    /* nA and nB are numbers less than the length of input A and B.             */
    /*     One may calculate the result not for all elements of input A and B.  */
                                                                 
    /* A[0:(nA-1)] is the result set A intersect B, which overwrites input A.                   */ 
    /*     (1) nonempty or empty                                                                */
    /*     (2) sorted                                                                           */
    /*     (3) no replcates for some elements.                                                  */    
    /*     So, if you want to keep the input A, do not use this function.                       */
    /* nA is the number of the result set A intersect B, which overwrites the length of input A.*/
    
    int lower = 0, upper, j = 0, i, midpoint;

    R_isort(A, nA); 
    R_isort(B, nB); 
    
    nA = nA-1;
    nB = nB-1;

    if (nA <= nB) {
        for (i=0; i<=nA; i++) 
        { 
            if (A[i] >= B[lower] && A[i] <= B[nB]) 
            {
                upper = nB;     
                if (upper == lower) 
                {
                    if (A[i] == B[lower]) 
                    {
                        A[j] = A[i];
                        j = j+1;
                        break;
                    }            
                } else if (upper == lower + 1) 
                {
                    if (A[i] == B[lower]) 
                    {
                        lower = lower + 1;
                        A[j] = A[i];
                        j = j+1;
                    } else if (A[i] == B[upper]) 
                    {
                        A[j] = A[i]; 
                        j = j+1;
                        break;
                    }           
                }
                while ((upper - lower) > 1)
                {
                    midpoint = (lower + upper) >> 1;
                    if (A[i] == B[lower]) 
                    {
                        lower = lower + 1;
                        A[j] = A[i]; 
                        j = j+1;
                        break;
                    } else if (A[i] == B[upper]) 
                    {
                        lower = upper;
                        A[j] = A[i];
                        j = j+1;
                        break;
                    } else if (A[i] == B[midpoint]) 
                    {
                        lower = midpoint + 1;
                        A[j] = A[i]; 
                        j = j+1;
                        break;
                    } else if (A[i] > B[midpoint])
                        lower = midpoint; 
                    else 
                        upper = midpoint; 
                }
            }
        }
        nA = j;
    } else 
    {
        for (i=0; i<=nB; i++) { 
            if (B[i] >= A[lower] && B[i] <= A[nA]) 
            {
                upper = nA;     
                if (upper == lower) 
                {
                    if (B[i] == A[lower]) {
                        A[j] = B[i];
                        j = j+1;
                        break;
                    }            
                } else if (upper == lower + 1) 
                {
                    if (B[i] == A[lower]) 
                    {
                        lower = lower + 1;
                        A[j] = B[i];
                        j = j+1;
                    } else if (B[i] == A[upper]) 
                    {
                        A[j] = B[i];
                        j = j+1;
                        break;
                    }           
                }
                while ((upper - lower) > 1)
                {
                    midpoint = (lower + upper) >> 1;
                    if (B[i] == A[lower]) 
                    {
                        lower = lower + 1;
                        A[j] = B[i];
                        j = j+1;
                        break;
                    } else if (B[i] == A[upper]) 
                    {
                        lower = upper;
                        A[j] = B[i];
                        j = j+1;
                        break;
                    } else if (B[i] == A[midpoint]) 
                    {
                        lower = midpoint + 1;
                        A[j] = B[i];
                        j = j+1;
                        break;
                    } else if (B[i] > A[midpoint])
                        lower = midpoint;
                    else 
                        upper = midpoint;
                }
            }
        }
        nA = j; 
    }
    return(nA);
}

void outerCC(int *A, int *nA, int *B, int *nB, int *OUTER, int *nOUTER) 
{
    *nOUTER = outerC(A, *nA, B, *nB, OUTER);
}   
static int outerC(int *A, int nA, int *B, int nB, int *OUTER)
{

    /* Finding A[i] + B[j] for all i, j.                                        */
    /* A and B are the set of integers which is                                 */
    /*     (1) nonempty                                                         */
    /*     (2) unsorted or sorted                                               */
    /*     (3) replcates for some elements.                                     */
    /* nA and nB are numbers less than the length of input A and B.             */
    /*     One may calculate the result not for all elements of input A and B.  */
                                                                 
    /* OUTER[0:(nOUTER-1)] is the sorted result {A[i] + B[j] for all i, j} - A. */ 
    /*     (1) nonempty                                                         */
    /*     (2) sorted                                                           */
    /*     (3) no replcates for some elements.                                  */
    /* nOUTER is the number of sorted result {A[i] + B[j] for all i, j} - A.    */ 
    
    int i, j, nOUTER;

    for (i=0; i<nA; i++)
        for (j=0; j<nB; j++)
            OUTER[i+nA*j] = A[i] + B[j];
        
    R_isort(OUTER, nA*nB);  
    nOUTER = uniqueC(OUTER, nA*nB);
    nOUTER = setdiffC(OUTER, nOUTER, A, nA);
    return(nOUTER);
}

void extrema2dVC(double *zz, int *nnrow, int *nncol,
              int *maxindex, int *nnmax, int *ttotalmax, int *minindex, int *nnmin, int *ttotalmin, int *vertex) 
{
    int nrow = *nnrow, ncol = *nncol, nvert = nrow * ncol;
    int i, j, nbound=2*(nrow+ncol-2), ninterior, npatch, nneighbor, nmax, nmin, totalmax, totalmin;
    
    /* int nb[8]    = {1, ncol-2,      1, nrow-2,           1,     ncol-2,              1,        nrow-2};  */
    int cumnb[8] = {0,      1, ncol-1,   ncol, ncol+nrow-2, ncol+nrow-1, 2*ncol+nrow-3, 2*ncol+nrow-2};    
 
    int *vertexbound = (int *) R_alloc(nbound, sizeof(int));
            
    vertexbound[cumnb[0]] = 0;   
    for (i=1; i<=(ncol-2); i++) vertexbound[cumnb[0]+i] = nrow * i; 
    vertexbound[cumnb[2]] = nrow * (ncol-1);   
    for (i=1; i<=(nrow-2); i++) vertexbound[cumnb[2]+i] = i + nrow * (ncol - 1); 
    vertexbound[cumnb[4]] = nrow - 1 + nrow * (ncol - 1);   
    for (i=1; i<=(ncol-2); i++) vertexbound[cumnb[4]+i] = nrow - 1 + nrow * i; 
    vertexbound[cumnb[6]] = nrow - 1;   
    for (i=1; i<=(nrow-2); i++) vertexbound[cumnb[6]+i] = i; 
                  
    /*int *vertex = (int *) R_alloc(nvert, sizeof(int)); */
    int *uvertex = (int *) R_alloc(nvert, sizeof(int));
    
    /*  if(vertex != NULL) Rprintf("%c memory for vertex \n", 1);   */
    /*  if(uvertex != NULL) Rprintf("%c memory for unique \n", 1);  */
    /*  Rprintf("%d\n", i); */
            
    /* Finding equivalence class induced by neighborhood relation */ 
    /* Neighborhood is defined if horizontal, vertical, diagonal adjacent pixels are the same. */
    /* Partition the image by Rem's algorithm */
     
    findclassC(zz, nrow, ncol, vertex);

    for (i=0; i<nvert; i++)
        uvertex[i] = vertex[i];
        
    /* Finding neighbors of patches including boundary vertex. */

    for (i=0; i<nbound; i++) 
        vertexbound[i] = vertex[vertexbound[i]];

    nbound = uniqueC(vertexbound, nbound);
    
    nmax = 0;
    nmin = 0;
    totalmax = 0;
    totalmin = 0;
   
    int *patch = (int *) R_alloc(nvert, sizeof(int));
    int *neighbor = (int *) R_alloc(nvert, sizeof(int));
       
    for (i=0; i<nbound; i++) {
        npatch = findindexC(vertex, nvert, vertexbound[i], patch);       
        nneighbor = neighborC(patch, npatch, nrow, ncol, 0, neighbor);
        if (nneighbor > 0) {
            if(allC(zz, neighbor, nneighbor, zz[patch[0]]) == nneighbor) {
                for (j=0; j<npatch; j++) 
                    maxindex[totalmax + j] = patch[j];
                    
                totalmax = totalmax + npatch;
                nmax = nmax + 1;
            }
            else if(allC(zz, neighbor, nneighbor, zz[patch[0]]) == -nneighbor) {
                for (j=0; j<npatch; j++) 
                    minindex[totalmin + j] = patch[j];
                                  
                totalmin = totalmin + npatch;
                nmin = nmin + 1;                
            }
        }  
    }


    ninterior = uniqueC(uvertex, nvert);
    ninterior = setdiffC(uvertex, ninterior, vertexbound, nbound);   
    if (ninterior != 0) {

        for (i=0; i<ninterior; i++) {
            npatch = findindexC(vertex, nvert, uvertex[i], patch);
            nneighbor = neighborC(patch, npatch, nrow, ncol, 1, neighbor);
            if (nneighbor > 0) {
                if(allC(zz, neighbor, nneighbor, zz[patch[0]]) == nneighbor) {
                    for (j=0; j<npatch; j++) 
                        maxindex[totalmax + j] = patch[j];
                        
                    totalmax = totalmax + npatch;
                    nmax = nmax + 1;
                }
                else if(allC(zz, neighbor, nneighbor, zz[patch[0]]) == -nneighbor) {
                    for (j=0; j<npatch; j++) 
                        minindex[totalmin + j] = patch[j];
                        
                    totalmin = totalmin + npatch;
                    nmin = nmin + 1;                
                }
            }
        }   
    }
   
    *nnmax = nmax;
    *nnmin = nmin;
    *ttotalmax = totalmax;
    *ttotalmin = totalmin;
}

#include <R_ext/Rdynload.h>

static const R_CMethodDef cMethods[] = {
    {"extrema2dC", (DL_FUNC) &extrema2dC, 9},
    {NULL, NULL, 0}
};

void R_init_EMD(DllInfo *info)
{
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}

