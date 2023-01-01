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

1. In each step calculate r_k = Ax_k − b
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

Tento repozitár obsahuje kód implementující Jacobi a Gauss-Seidel metodu pro vyřešení konkrétního typu lineárních rovnic.

**Obsah:**
- [Instrukce](#instrukce)
- [Iterativní metoda](#iterativni-metoda)
  - [Jacobi metoda](#jacobi-metoda)
  - [Gauss-Seidel metoda](#gauss-seidel-metoda)
- [Záruky konvergence](#záruky-konvergence)
- [Výsledky](#výsledky)

# Instrukce

Zde je pomocná zpráva skriptu:

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

Talže, pokud chcete použít Gausse-Seidel metodu s gamma 2 a iterativním zpřesněním:

```
python3 lineqs.py -g 2 -m gauss_seidel -r true
```

Nebo Jacobi metodu s gamma 10 bez zpřesnění:

```
python3 lineqs.py -g 10 -m jacobi
```
# Iterativní metoda

Cílem je sestrojit algoritmus, který řeší Ax = b. Tento přístup volí x_0 jako nulový vektor a snaží se iterativně změnit tento vektor, tak aby byl čím dál blíž x.

Pro to aby jsme toho dosáhli zvolme nejdříve náhodnou regulární matici Q (Jacobi a Gauss-Seidel metoda se liší pouze s výběrem této matice Q). Poté pro spočítání x_k:

x_k = Q^-1((Q - A)x_k-1 - b)

Pokud posloupnost x_k konverguje, máme řešení.

Těžká čést tohoto přístupu je jak najít Q. Pro pomoc s výběrem lze upravit rovnici nahoře takto:

x_k - x = (E - Q^-1*A)(x_k-1 - x) 

Kde e_k = x_k - x se nazývá vektorem chyby a E - Q^-1*A = W.

Lze si všimnout, že  e_k = W^1 * e_k-1 = W^2 * e_k-2 = ... = W^k * e_0

Z čehož plyne, že lim k-> inf W = 0, což je ekvivalentní s p(W) < 1, kde p(M je spektrální poloměr matice M.

Dále můžeme udělat konvergenci rychlejší volbout W takovou, že ||W|| je co nejmenší (||M|| je maticová norma).

Pro shrnutí, chceme aby matice W (W = E - Q^-1*A) měla spektrální poloměr < 1 a aby její norma byla co nejmenší. Jacobi a Gauss-Seidel metoda použili tyto "pravidla" pro konstrukci jejich Q matic.

## Jacobi metoda

Jacobi metoda nastavuje Q = D, kde D je nulová matice s diagonálou a.

## Gauss-Seidel metoda

Gauss-Seidel metoda nastavuje Q = D + L, kde D je nulová matice s diagonálou A a L je nižší trojúhelník matice A - D.

# Záruky konvergence

Ve skriptu jsou čtyři kontroly, které mohou zaručit teoretickou konvergence x_k - dvě kontroly pro každou metodu.

První dvě kontrolují jestli je spektrální poloměr W menší než 1. Pokud je, zaručuje to teoretickou konvergenci.

Poté je v Jacobi metodě zkontrolováno zda je A ostře diagonálně dominantní. Pokud je, zaručuje to teoretickou konvergenci.

Nakonec, v Gauss-Seidel metodě je zkontrolováno, zdali je A symetrická, pozitivně deifinitivní a má pozitivní diagonálu. Pokud ano, zaručuje to konvergenci.

# Výsledky

Zde jsou dvě tabulky, první používá standardní iteratní metodu vysvětlenou nahoře, ta druhá používá iterativní zpřesnění:

1. V každém kroku vypočítat r_k = Ax_k − b
2. Vyřešit Ay_k = r_k
3. Změnit x_k : x_k = x_k - y_k 

Standardní iterativní metoda (NaN znamená že metoda nekonvergovala):

|Metoda      |Gamma|#Kroky|Error                 |
|:-----------|:---:|:----:|:--------------------:|
|Jacobi      |10   |8     |2.1601349385807865e-06|
|Gauss-Seidel|10   |6     |1.5714747643696884e-06|
|Jacobi      |2    |782   |9.968534524306131e-06 |
|Gauss-Seidel|2    |393   |9.807646243637133e-06 |
|Jacobi      |4/5  | Nan  | Nan                  |
|Gauss-Seidel|4/5  | Nan  | Nan                  |

Iterativní zpřesnění:

|Metoda      |Gamma|#Kroky|Error                  |
|:-----------|:---:|:----:|:---------------------:|
|Jacobi      |10   | 1    |0                      |
|Gauss-Seidel|10   | 1    |8.487761295006218e-17  |
|Jacobi      |2    | 1    |4.839349969133126e-16  |
|Gauss-Seidel|2    | 1    |4.1540741810552243e-16 |
|Jacobi      |4/5  | 1    |5.1201905234891505e-16 |
|Gauss-Seidel|4/5  | 1    |3.659787256322971e-15  |
