/*
Home
Discussions
Write at Opengenus IQ
HOME
😊 JOIN OUR INTERNSHIP 🎓
MACHINE LEARNING 🤖
ALGORITHMS
DATA STRUCTURES
LEARN C++
MARKDOWN GUIDE
Toom-Cook method for multiplication
Join the strongest 💪 computer science community in the World for free.
multiplication integer multiplication toom cook algorithm algorithm
 Feature Toggle Rollout features continuously. Why wait for a new build? Flip a toggle.ethical ad by CodeFund 
Reading time: 25 minutes | Coding time: 10 minutes

Multiplication of two n-digits integers has time complexity at worst O(n^2).Toom-Cook algorithm is an algorithm for multiplying two n digit numbers in Θ(c(k)n^e) time complexity , where e = log(2k − 1) / log(k), n^e is the time spent on sub-multiplications, and c is the time spent on additions and multiplication by small constants.

The idea is based on divide-and-conquer technique. Given two large integers, a and b, Toom–Cook splits up a and b into k smaller parts each of length l, and performs operations on the parts. As k grows, one may combine many of the multiplication sub-operations, thus reducing the overall complexity of the algorithm. The multiplication sub-operations can then be computed recursively using Toom–Cook multiplication again, and so on.

Toom-Cook 3-Way Multiplication
Explanation
Toom cook algorithm is the advanced approach for splitting the numbers into parts. Toom cook n way reduces the product to 2*(n)-1 multiplications.Where n stands for 3.Let the operands considered are split into 3 pieces of equal length. The parts are written in terms of polynomials:

X(t) = (X2)t^2 + X1(t) + X0
Y(t) = (Y2)t^2 + Y1(t) + Y0

base B = b is choosen , such that the number of digits of both x and m in base B is at most k.Here,we take Base B = b^2 = 108.Then we will seperate x and y into base B digits x,y.These 2 equations are multiplied to form w(t)=x(t) * y(t).

W(t) = w4 ∗ t4 + w3 ∗ t3 + w2 ∗ t2 + w1 ∗ t + w0

The final w(t) is calculated through the value of t, although the final step is going to be the addition.X(t) and Y(t) are calculated and multiplied by choosing some set of points, forming w(t).Let the following points be (0,1,2,-1,inf)

t=0 x0 * y0
t=1 (x2+x1+x0) * (y2+y1+y0)
t=-1 (x2-x1+x0) * (y2-y1+y0)
t = 2(4 ∗ x2 + 2 ∗ x1 + x0) ∗ (4 ∗ y2 + 2 ∗ y1 + y0)
t=inf x2 * y2
Then, the value of those combinations is calculated through

W(0) = w0
W(1) = w4 + w3 + w2 + w1 + w0
W(-1) = w4 - w3 + w2 - w1 + w0
W(2) = 16 * w4 + 8 * w3 + 4 * w2 + 2 * w1 + w0
W(inf) = w4
Toom-3 running time is significally O(n^1.465).

Pseudocode
1.Input : TwointegersAandBaregivenwhere0 < A, B, < Xn

2.Output : AB = c0 + c1x^k + c2^x2k + c3x^3k + c4^x4k whenk = n/(k = numberofsplits)

3.Here A = x0 + x1t + x2t^2, B = y0 + y1t + y2t^2 here t = X^k

t=(0,1,2,-1,inf)

t=0 x0 * y0, which gives w0 immediately

t=1 (x2+x1+x0) * (y2+y1+y0)

t=-1 (x2-x1+x0) * (y2-y1+y0)

Example
Consider 2 numbers 831275469 , 897512436.Both the numbers are split into 3 lim bs each of length 3 digits written into polynomials.

P(x) = a2 * x^2 + a1 * x + a0(a2 = 831; a1 = 275; a0 = 469)
Q(x) = b2 * x^2 + b1 * x + b0(b2 = 897; b1 = 512; b0 = 436)

p(x) = 831x^2 + 275x + 469
q(x) = 897x^2 + 512x + 436

We write, p(x)q(x) = r(x).
(831x^2 +275x+469)(897x^2 +512x+436) = ax^4 +bx^3 +cx^2+dx+e= r(x).

Then we substitute the values of x for set of equations.Lets the set of points to be substituted be -2,-1,0,1,2.

For x=0:
p(0)q(0) = r(0).By solving this equation we get e = 204484.
For x=1:
p(1)q(1) = r(1).By solving this equation we get a+b+c+d+e = 2905875 ||eqation 1.
For x=-1:
p(-1)q(-1) = r(-1).By solving this equation we get a-b+c-d+e = 841525 ||equation 2.
For x=2:
p(2)q(2) = r(2).By solving this equation we get 16a+8b+4c+2d+e= 21923464 ||equation 3.
For x=-2:
p(-2)q(-2) = r(-2).By solving this equation we get 16a-8b+4c-2d+e= 1027116 ||equation 4.
The equations we got by those points are:

e = 204484
a+ b+ c+ d+e = 2905875
a-b+ c- d+e = 841525
16a+8b+4c+2d+e = 21923464
16a-8b+4c-2d+e = 1027116
By solving the above equation we get
a = 382828 , c = -1286387 , b = 1397304 , d = -365129 , e=204484.
The Final product is obtained by adding the numbers by shifting,The result is:38296772326152333484

Complexity
We see that the k-way split requires 2 k - 1 multiplications.This means that instead of k2 multiplications we do only 2 k - 1.The recurrence equation for the total number of multiplication is given by:T(n) = (2 k - 1)T(n/k) + O(n).

Time Complexity of Toom–Cook multiplication will be:

T(n) = (2 k - 1)^logk(n) = n logk (2 k-1)

Implementations
C++
*/

//Toom-Cook 3-Way Multiplication
#include <cstdlib>   // for rand()
#include <iostream>  // for cout
#include <math.h>    // for pow()
#include <stdio.h>   // for printf()
//#define TEST      
#define D_MAX 729  
#define D     729  
using namespace std;
class Calc
{
    int A[D];  
    int B[D];  
#ifdef TEST
    int cnt_mul;       
    clock_t t1, t2;    
    double tt;         
#endif
    public:
        Calc();                                         
        void calcToomCook();                                
    private:
        void multiplyNormal(int *, int *, int, int *);      
        void multiplyToomCook3(int *, int *, int , int *);  
        void doCarry(int *, int);                           
        void display(int *, int *, int *);          
};
Calc::Calc()
{
     int i;  
    for (i = 0; i < D; i++) 
    {
        A[i] = rand() % 10;
        B[i] = rand() % 10;
    }
}
void Calc::calcToomCook()
{
    int a[D_MAX];       
    int b[D_MAX];       
    int z[D_MAX * 2];  
    int i;              
#ifdef TEST
    t1 = clock();  
    for (int l = 0; l < 1000; l++) {
        cnt_mul = 0;  
#endif
        for (i = 0; i < D; i++)
    {
        a[i] = A[i];
        b[i] = B[i];
    }
        for (i = D; i < D_MAX; i++)
    {
        a[i] = 0;
        b[i] = 0;
    }
    multiplyToomCook3(a, b, D_MAX, z);
    doCarry(z, D_MAX * 2);
#ifdef TEST
    }
    t2 = clock();  
    tt = (double)(t2 - t1) / CLOCKS_PER_SEC; 
#endif
display(a, b, z);
}
/*
 * Multiplication Standard
 */
void Calc::multiplyNormal(int *a, int *b, int tLen, int *z)
{
    int i, j;  
   for(i = 0; i < tLen * 2; i++) z[i] = 0;
   for (j = 0; j < tLen; j++) {
        for (i = 0; i < tLen; i++) {
            z[j + i] += a[i] * b[j];
#ifdef TEST
            cnt_mul++; 
#endif
        }
    }
}
void Calc::multiplyToomCook3(int *a, int *b, int tLen, int *z)
{
    int *a0 = &a[0];                // Multiplicand / right side array pointer
    int *a1 = &a[tLen / 3];         // Multiplicand / central array pointer
    int *a2 = &a[tLen * 2/ 3];      // Multiplicand / left side array pointer
    int *b0 = &b[0];                // Multiplier / right side array pointer
    int *b1 = &b[tLen / 3];         // Multiplier / central array pointer
    int *b2 = &b[tLen * 2/ 3];      // Multiplier / left side array pointer
    int *c0 = &z[(tLen / 3) *  0];  // c0     
    int *c2 = &z[(tLen / 3) *  2];  // c2     
    int *c4 = &z[(tLen / 3) *  4];  // c4    
    int c1      [(tLen / 3) * 2];   // c1    
    int c3      [(tLen / 3) * 2];   // c3     
    int a_m2    [tLen / 3];         // a( -2) 
    int a_m1    [tLen / 3];         // a( -1) 
    int a_0     [tLen / 3];         // a(  0) 
    int a_1     [tLen / 3];         // a(  1) 
    int a_inf   [tLen / 3];         // a(inf) 
    int b_m2    [tLen / 3];         // b( -2) 
    int b_m1    [tLen / 3];         // b( -1) 
    int b_0     [tLen / 3];         // b(  0) 
    int b_1     [tLen / 3];         // b(  1)
    int b_inf   [tLen / 3];         // b(inf) 
    int c_m2    [(tLen / 3) * 2];   // c( -2)
    int c_m1    [(tLen / 3) * 2];   // c( -1)
    int c_0     [(tLen / 3) * 2];   // c(  0) 
    int c_1     [(tLen / 3) * 2];   // c(  1) 
    int c_inf   [(tLen / 3) * 2];   // c(inf) 
    int i;                          
    if (tLen <= 9) 
    {
        multiplyNormal(a, b, tLen, z);
        return;
    }
    // ==== a(-2) = 4 * a2 - 2 * a1 + a0, b(1) = 4 * b2 - 2 * b1 + b0
    for(i = 0; i < tLen / 3; i++)
    {
        a_m2[i] = (a2[i] << 2) - (a1[i] << 1) + a0[i];
        b_m2[i] = (b2[i] << 2) - (b1[i] << 1) + b0[i];
    }
    // ==== a(-1) = a2 - a1 + a0, b(1) = b2 - b1 + b0
    for(i = 0; i < tLen / 3; i++) 
    {
        a_m1[i] = a2[i] - a1[i] + a0[i];
        b_m1[i] = b2[i] - b1[i] + b0[i];
    }
    // ==== a(0) = a0, b(0) = b0
    for(i = 0; i < tLen / 3; i++) {
        a_0[i] = a0[i];
        b_0[i] = b0[i];
    }
    // ==== a(1) = a2 + a1 + a0, b(1) = b2 + b1 + b0
    for(i = 0; i < tLen / 3; i++) {
        a_1[i] = a2[i] + a1[i] + a0[i];
        b_1[i] = b2[i] + b1[i] + b0[i];
    }
    // ==== a(inf) = a2, b(inf) = b2
    for(i = 0; i < tLen / 3; i++) {
        a_inf[i] = a2[i];
        b_inf[i] = b2[i];
    }
    // ==== c(-2) = a(-2) * b(-2)
    multiplyToomCook3(a_m2,  b_m2,  tLen / 3, c_m2 );
    // ==== c(-1) = a(-1) * b(-1)
    multiplyToomCook3(a_m1,  b_m1,  tLen / 3, c_m1 );
    // ==== c(0) = a(0) * b(0)
    multiplyToomCook3(a_0,   b_0,   tLen / 3, c_0  );
    // ==== c(1) = a(1) * b(1)
    multiplyToomCook3(a_1,   b_1,   tLen / 3, c_1  );
    // ==== c(inf) = a(inf) * b(inf)
    multiplyToomCook3(a_inf, b_inf, tLen / 3, c_inf);
    // ==== c4 = 6 * c(inf) / 6
    for(i = 0; i < (tLen / 3) * 2; i++)
        c4[i] = c_inf[i];
        // ==== c3 = -c(-2) + 3 * c(-1) - 3 * c(0) + c(1) + 12 * c(inf) / 6
    for(i = 0; i < (tLen / 3) * 2; i++)
    {
        c3[i]  = -c_m2[i];
        c3[i] += (c_m1[i] << 1) + c_m1[i];
        c3[i] -= (c_0[i] << 1) + c_0[i];
        c3[i] += c_1[i];
        c3[i] += (c_inf[i] << 3) + (c_inf[i] << 2);
        c3[i] /= 6;
    }
    // ==== c2 = 3 * c(-1) - 6 * c(0) + 3 * c(1) - 6 * c(inf) / 6
    for(i = 0; i < (tLen / 3) * 2; i++) {
        c2[i]  = (c_m1[i] << 1) + c_m1[i];
        c2[i] -= (c_0[i] << 2) + (c_0[i] << 1);
        c2[i] += (c_1[i] << 1) + c_1[i];
        c2[i] -= (c_inf[i] << 2) + (c_inf[i] << 1);
        c2[i] /= 6;
    }
    // ==== c1 = c(-2) - 6 * c(-1) + 3 * c(0) + 2 * c(1) - 12 * c(inf) / 6
    for(i = 0; i < (tLen / 3) * 2; i++) {
        c1[i]  = c_m2[i];
        c1[i] -= (c_m1[i] << 2) + (c_m1[i] << 1);
        c1[i] += (c_0[i] << 1) + c_0[i];
        c1[i] += (c_1[i] << 1);
        c1[i] -= (c_inf[i] << 3) + (c_inf[i] << 2);
        c1[i] /= 6;
    }
    // ==== c0 = 6 * c(0) / 6
    for(i = 0; i < (tLen / 3) * 2; i++)
        c0[i] = c_0[i];
    // ==== z = c4 * x^4 + c3 * x^3 + c2 * x^2 + c1 * x + c0
    for(i = 0; i < (tLen / 3) * 2; i++) z[i + tLen / 3] += c1[i];
    for(i = 0; i < (tLen / 3) * 2; i++) z[i + (tLen / 3) * 3] += c3[i];
}
void Calc::doCarry(int *a, int tLen) {
    int cr;  
    int i;   
    cr = 0;
    for(i = 0; i < tLen; i++) {
        a[i] += cr;
        if(a[i] < 0) {
            cr = -(-(a[i] + 1) / 10 + 1);
        } else {
            cr = a[i] / 10;
        }
        a[i] -= cr * 10;
    }
    // Overflow
    if (cr != 0) printf("[ OVERFLOW!! ] %d\n", cr);
}
/*
 * Result output
 */
void Calc::display(int *a, int *b, int *z)
{
    int i; 
    int aLen = D_MAX, bLen = D_MAX, zLen = D_MAX * 2;
    while (a[aLen - 1] == 0) if (a[aLen - 1] == 0) aLen--;
    while (b[bLen - 1] == 0) if (b[bLen - 1] == 0) bLen--;
    while (z[zLen - 1] == 0) if (z[zLen - 1] == 0) zLen--;
    // a 
    printf("a =\n");
    for (i = aLen - 1; i >= 0; i--) {
        printf("%d", a[i]);
        if ((aLen - i) % 10 == 0) printf(" ");
        if ((aLen - i) % 50 == 0) printf("\n");
    }
    printf("\n");
    // b 
    printf("b =\n");
    for (i = bLen - 1; i >= 0; i--) {
        printf("%d", b[i]);
        if ((bLen - i) % 10 == 0) printf(" ");
        if ((bLen - i) % 50 == 0) printf("\n");
    }
    printf("\n");
    // z 
    printf("z =\n");
    for (i = zLen - 1; i >= 0; i--) {
        printf("%d", z[i]);
        if ((zLen - i) % 10 == 0) printf(" ");
        if ((zLen - i) % 50 == 0) printf("\n");
    }
    printf("\n\n");
#ifdef TEST
    printf("Counts of multiply / 1 loop = %d\n", cnt_mul);     // Multiplication count
    printf("Total time of all loops     = %f seconds\n", tt);  // processing time
#endif
}
int main()
{
    try
    {
        Calc objCalc;
        objCalc.calcToomCook();
    }
    catch (...) {
        cout << "Exception occurred!" << endl;
        return -1;
    }
    return 0;
}

/*
C++Copy
Applications
There are various areas where this algorithm application is done, which involves multiplication of large integers.

This method is used in McEliece Cryptosystems to overcome the drawbacks in terms of size of the encrypted key and transmission rate.

Big Number arithmetic.

In Cryptographic algorithms especially for reducing the complexity in encoding and decoding of the keys like ElGamal, RSA, Elliptical Curvecryptosystems,and Diffie Hellman key exchange protocol.

Calculation of mathematical constants like Pi, e etc.

Kyatham Srikanth
Read more posts by this author.

Read More
OpenGenus Foundation OpenGenus Foundation
 Tags
multiplication
integer multiplication
toom cook algorithm
algorithm

— OpenGenus IQ: Learn Computer Science —
multiplication
Karatsuba Algorithm (for fast integer multiplication)
1 post →
C
Semaphore in C
Semaphore is a data handling technique which is very useful in process synchronization and multithreading. We used the POSIX semaphore library and use the functions sem_wait, sem_post, sem_open and sem_init for implementing semaphore in C.

OpenGenus Foundation OPENGENUS FOUNDATION
MATRIX MULTIPLICATION
Cannon’s algorithm for distributed matrix multiplication
Cannon's algorithm is a distributed algorithm for matrix multiplication for two-dimensional meshes. It is especially suitable for computers laid out in an N × N mesh. The main advantage of the algorithm is that its storage requirements remain constant and are independent of the number of processors.

KYATHAM SRIKANTH
OpenGenus IQ: Learn Computer Science icon
OpenGenus IQ: Learn Computer Science
—
Toom-Cook method for multiplication
Share this


OpenGenus IQ: Learn Computer Science © 2019 All rights reserved ™
Top Posts
Facebook
Twitter
*/
