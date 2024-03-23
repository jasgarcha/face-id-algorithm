# face-id-algorithm
A work-in-progress implementation of an algorithm for biometric facial identification/recognition. Currently, it is a little less than an image comparator.

## Logic
Algorithm:
Start with 2 matrices represented by 2-dimensional arrays.
1) Matrices -> vectors (column vector size = matrix number of rows * matrix number of columns).

2) Center, normalize, linearize matrices.

Let `S` be a set of `M` face images. `S = {T1, T2,...TM}`.

2) Compute the mean image: mean image = `(1/M)*T1+(1/M)*T2+...+(1/M)*TM`.

3) Compute the difference between the input image and the mean image: `T1-mean image, T2-mean image,...,Ti-mean image`.


`A = [Image1AndMeanImageDifferenceImage2AndMeanImageDifference]`

4) Work around: Find the eigenvalues and eigenvectors of `A` transpose times `A`.

 4.1. Get eigenvalues and eigenvectors of `M`, not the work around `ATA`. `u = Av`.

5) Compute normalized eigenvector.

6) New representation: Representing a face onto the basis of the normalized eigenvector.

`I = w*u` where `w = normalizeUTranspose*ImageAndMeanImageDifference`
...

## Todo
Implement rest of algorithm:

Add scalars that acts as a weight `w_j` to the eigenvector in the linear combination expression for representing the original image.

7) Feature vector = set of weight = `{w1, w2,..., wk}`.

For every image vector, compare vectors.

Compare against a third image. 

Compare `wT*InputImageAndMeanImageDifference3` against `wT*InputImageAndMeanImageDifference1` and `wT*InputImageAndMeanImageDifference2`.

Add decision function.

## Build
The third-party Apache Commons Mathematics library is required for linear algebra operations.