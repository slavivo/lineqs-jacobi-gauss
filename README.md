# English

# Jacobi and Gauss-Seidel methods

This repository contains code that implements Jacobi and Gauss-Seidel methods to solve a certain type of linear equations.

**Table of contents:**
- [Instructions](#instructions)
- [Iterative approach](#iterative-approach)
  - [Jacobi method](#jacobi-method)
  - [Gauss-Seidel method](#gauss-seidel-method)
- [Guarantees of convergence](#guarantees-of-convergence)
- [Results](#results)

# Instructions

Here is the help message of the script:

```
python3 lineqs.py --help
usage: lineqs.py [-h] [-g GAMMA] [-m {jacobi,gauss_seidel}] [-r REFINEMENT]

optional arguments:
  -h, --help            show this help message and exit
  -g GAMMA, --gamma GAMMA
                        Parameter which is used in construction of matrix A
                        and vector b
  -m {jacobi,gauss_seidel}, --method {jacobi,gauss_seidel}
                        Method to solve the system of linear equations
  -r REFINEMENT, --refinement REFINEMENT
                        Use iterative refinement
```

So, if you wanted to use the Gauss-Seidel method with gamma of 2 with refinement:

```
python3 lineqs.py -g 2 -m gauss_seidel -r true
```

Or the Jacobi method with gamma 10 without refinement:

```
python3 lineqs.py -g 10 -m jacobi
```
# Iterative approach

The goal is to construct an algorithm that solves Ax = b. This approach chooses an x_0 as a null vector and seeks to iteratively change it to be closer and closer to x.

To do this, let's first choose a random regular matrix Q (the Jacobi and Gauss-Seidel methods have only one difference and that's how they choose this Q). Then to calculate x_k :

x_k = Q^-1((Q - A)x_k-1 - b)

If the x_k sequence converges we have a solution.

The difficult part of this approach is how to choose Q. To help with that we can modify the above equation to this:

x_k - x = (E - Q^-1*A)(x_k-1 - x) 

Where e_k = x_k - x is called the vector error and E - Q^-1*A = W.

We can see that: e_k = W^1 * e_k-1 = W^2 * e_k-2 = ... = W^k * e_0

So that means that lim k-> inf W = 0, which is equal to: p(W) < 1, where p(M) is spectral radius of matrix M.

We can also make the convergence faster by choosing W so that ||W|| is as small as possible (||M|| is the matrix norm).

To summarise, we want the W matrix (W = E - Q^-1*A) to have spectral redius < 1 and for ||W|| to be as small as possible. The Jacobi and Gauss-Seidel methods used these "rules" to construct their Q matrices.

## Jacobi method

The Jacobi method sets Q = D, where D is a zero filled matrix with diagonal of A.

## Gauss-Seidel method

The Gauss-Seidel method sets Q = D + L, where D is a zero filled matrix with diagonal of A and L is the lower triangle of A - D.

# Guarantees of convergence

In the script there are four checks that can guarantee theoretical convergence of x_k - two for each method.

The first two check if the spectral radius of W is lesser than 1. If it is, it guaratees theoretical convergence.

Then in the Jacobi we check if A is is sharply diagonally dominant. If it is, it guarantees theoretical convergence.

Lastly, in Gauss-Seidel method we check if A is symmetric and positive definite and has positive diagonal. If it is, it guarantees theoretical convergence.

# Results

Here are two tables, the first one uses the above described iterative approach, but the other uses so called iterative refinement:

1. In each step calculate r_k = Ax_k ??? b
2. Solve Ay_k = r_k
3. Change x_k : x_k = x_k - y_k 

Default iterative approach (NaN means that the method did not converge):

|Method      |Gamma|#Steps|Error                 |
|:-----------|:---:|:----:|:--------------------:|
|Jacobi      |10   |8     |2.1601349385807865e-06|
|Gauss-Seidel|10   |6     |1.5714747643696884e-06|
|Jacobi      |2    |782   |9.968534524306131e-06 |
|Gauss-Seidel|2    |393   |9.807646243637133e-06 |
|Jacobi      |4/5  | Nan  | Nan                  |
|Gauss-Seidel|4/5  | Nan  | Nan                  |

Iterative refinement:

|Method      |Gamma|#Steps|Error                  |
|:-----------|:---:|:----:|:---------------------:|
|Jacobi      |10   | 1    |0                      |
|Gauss-Seidel|10   | 1    |8.487761295006218e-17  |
|Jacobi      |2    | 1    |4.839349969133126e-16  |
|Gauss-Seidel|2    | 1    |4.1540741810552243e-16 |
|Jacobi      |4/5  | 1    |5.1201905234891505e-16 |
|Gauss-Seidel|4/5  | 1    |3.659787256322971e-15  |

# Czech

# Jacobi a Gauss-Seidel methody

Tento repozit??r obsahuje k??d implementuj??c?? Jacobi a Gauss-Seidel metodu pro vy??e??en?? konkr??tn??ho typu line??rn??ch rovnic.

**Obsah:**
- [Instrukce](#instrukce)
- [Iterativn?? metoda](#iterativni-metoda)
  - [Jacobi metoda](#jacobi-metoda)
  - [Gauss-Seidel metoda](#gauss-seidel-metoda)
- [Z??ruky konvergence](#z??ruky-konvergence)
- [V??sledky](#v??sledky)

# Instrukce

Zde je pomocn?? zpr??va skriptu:

```
python3 lineqs.py --help
usage: lineqs.py [-h] [-g GAMMA] [-m {jacobi,gauss_seidel}] [-r REFINEMENT]

optional arguments:
  -h, --help            show this help message and exit
  -g GAMMA, --gamma GAMMA
                        Parameter which is used in construction of matrix A
                        and vector b
  -m {jacobi,gauss_seidel}, --method {jacobi,gauss_seidel}
                        Method to solve the system of linear equations
  -r REFINEMENT, --refinement REFINEMENT
                        Use iterative refinement
```

Tal??e, pokud chcete pou????t Gausse-Seidel metodu s gamma 2 a iterativn??m zp??esn??n??m:

```
python3 lineqs.py -g 2 -m gauss_seidel -r true
```

Nebo Jacobi metodu s gamma 10 bez zp??esn??n??:

```
python3 lineqs.py -g 10 -m jacobi
```
# Iterativn?? metoda

C??lem je sestrojit algoritmus, kter?? ??e???? Ax = b. Tento p????stup vol?? x_0 jako nulov?? vektor a sna???? se iterativn?? zm??nit tento vektor, tak aby byl ????m d??l bl???? x.

Pro to aby jsme toho dos??hli zvolme nejd????ve n??hodnou regul??rn?? matici Q (Jacobi a Gauss-Seidel metoda se li???? pouze s v??b??rem t??to matice Q). Pot?? pro spo????t??n?? x_k:

x_k = Q^-1((Q - A)x_k-1 - b)

Pokud posloupnost x_k konverguje, m??me ??e??en??.

T????k?? ????st tohoto p????stupu je jak naj??t Q. Pro pomoc s v??b??rem lze upravit rovnici naho??e takto:

x_k - x = (E - Q^-1*A)(x_k-1 - x) 

Kde e_k = x_k - x se naz??v?? vektorem chyby a E - Q^-1*A = W.

Lze si v??imnout, ??e  e_k = W^1 * e_k-1 = W^2 * e_k-2 = ... = W^k * e_0

Z ??eho?? plyne, ??e lim k-> inf W = 0, co?? je ekvivalentn?? s p(W) < 1, kde p(M je spektr??ln?? polom??r matice M.

D??le m????eme ud??lat konvergenci rychlej???? volbout W takovou, ??e ||W|| je co nejmen???? (||M|| je maticov?? norma).

Pro shrnut??, chceme aby matice W (W = E - Q^-1*A) m??la spektr??ln?? polom??r < 1 a aby jej?? norma byla co nejmen????. Jacobi a Gauss-Seidel metoda pou??ili tyto "pravidla" pro konstrukci jejich Q matic.

## Jacobi metoda

Jacobi metoda nastavuje Q = D, kde D je nulov?? matice s diagon??lou a.

## Gauss-Seidel metoda

Gauss-Seidel metoda nastavuje Q = D + L, kde D je nulov?? matice s diagon??lou A a L je ni?????? troj??heln??k matice A - D.

# Z??ruky konvergence

Ve skriptu jsou ??ty??i kontroly, kter?? mohou zaru??it teoretickou konvergence x_k - dv?? kontroly pro ka??dou metodu.

Prvn?? dv?? kontroluj?? jestli je spektr??ln?? polom??r W men???? ne?? 1. Pokud je, zaru??uje to teoretickou konvergenci.

Pot?? je v Jacobi metod?? zkontrolov??no zda je A ost??e diagon??ln?? dominantn??. Pokud je, zaru??uje to teoretickou konvergenci.

Nakonec, v Gauss-Seidel metod?? je zkontrolov??no, zdali je A symetrick??, pozitivn?? deifinitivn?? a m?? pozitivn?? diagon??lu. Pokud ano, zaru??uje to konvergenci.

# V??sledky

Zde jsou dv?? tabulky, prvn?? pou????v?? standardn?? iteratn?? metodu vysv??tlenou naho??e, ta druh?? pou????v?? iterativn?? zp??esn??n??:

1. V ka??d??m kroku vypo????tat r_k = Ax_k ??? b
2. Vy??e??it Ay_k = r_k
3. Zm??nit x_k : x_k = x_k - y_k 

Standardn?? iterativn?? metoda (NaN znamen?? ??e metoda nekonvergovala):

|Metoda      |Gamma|#Kroky|Error                 |
|:-----------|:---:|:----:|:--------------------:|
|Jacobi      |10   |8     |2.1601349385807865e-06|
|Gauss-Seidel|10   |6     |1.5714747643696884e-06|
|Jacobi      |2    |782   |9.968534524306131e-06 |
|Gauss-Seidel|2    |393   |9.807646243637133e-06 |
|Jacobi      |4/5  | Nan  | Nan                  |
|Gauss-Seidel|4/5  | Nan  | Nan                  |

Iterativn?? zp??esn??n??:

|Metoda      |Gamma|#Kroky|Error                  |
|:-----------|:---:|:----:|:---------------------:|
|Jacobi      |10   | 1    |0                      |
|Gauss-Seidel|10   | 1    |8.487761295006218e-17  |
|Jacobi      |2    | 1    |4.839349969133126e-16  |
|Gauss-Seidel|2    | 1    |4.1540741810552243e-16 |
|Jacobi      |4/5  | 1    |5.1201905234891505e-16 |
|Gauss-Seidel|4/5  | 1    |3.659787256322971e-15  |
