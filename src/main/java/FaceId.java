import org.apache.commons.math3.linear.*;

public class FaceId {
	public static void main(String[] args) {
		//double[][] image1 = {{1, 2}, {3, 4}};
		//double[][] image2 = {{5, 6}, {7, 8}};		
		double[][] image1 = {{5, 2}, {7, 9}};
		double[][] image2 = {{12, 7}, {2, 6}};		
		
		//Matrix -> Vector. Column vector size = Matrix number of rows * Matrix number of columns.
		
		//Center, normalize, linearize.
		double[][] T1 = new double[image1.length*image1[0].length][1]; //Square matrix, so image[0].length = image[i].length. 
		double[][] T2 = new double[image2.length*image2[0].length][1];		
		
		for(int i = 0; i < image1.length; i++)
			for(int j = 0; j < image1[i].length; j++)
				T1[i*image1.length+j][0] = image1[i][j];
		
		for(int i = 0; i < image2.length; i++)
			for(int j = 0; j < image2[i].length; j++)
				T2[i*image2.length+j][0] = image2[i][j];
		
		//double[][] T1 = {{1}, {2}, {3}, {4}}; Nx1 vector. 
		//double[][] T2 = {{5}, {6}, {7}, {8}}; Nx1 vector.	
		
		RealMatrix M1 = new Array2DRowRealMatrix(T1);
		RealMatrix M2 = new Array2DRowRealMatrix(T2);
		
		System.out.println("T1 = "+M1);
		System.out.println("T2 = "+M2);
		
		//Let S be a set of M face images. S = {T1, T2,...TM}
		RealMatrix[] S = {M1, M2};
		int M = S.length;
		
		//2. Compute the mean image: mean image = (1/M)*T1+(1/M)*T2+...+(1/M)*TM.
		RealMatrix meanImage = new Array2DRowRealMatrix(S[0].getRowDimension(), S[0].getColumnDimension());	
		for(int i = 0; i < M; i++) 
			meanImage = meanImage.add(S[i]); 
		meanImage = meanImage.scalarMultiply((double)1/M);
		
		//System.out.println("Mean Image = 1/"+M+"*("+M1+"+"+M2+")");
		System.out.println("Mean Image = "+meanImage);
		
		//3. Compute the difference between the input image and the mean image: T1-mean image, T2-mean image,...,Ti-mean image.
		RealMatrix Image1AndMeanImageDifference = M1.subtract(meanImage);
		RealMatrix Image2AndMeanImageDifference = M2.subtract(meanImage);
		
		System.out.println("Image1AndMeanImageDifference = "+Image1AndMeanImageDifference);
		System.out.println("Image2AndMeanImageDifference = "+Image2AndMeanImageDifference);
		
		//A = [Image1AndMeanImageDifferenceImage2AndMeanImageDifference]
		RealMatrix A = new Array2DRowRealMatrix(Image1AndMeanImageDifference.getRowDimension(), Image1AndMeanImageDifference.getColumnDimension()*2);
		for(int i = 0; i < Image1AndMeanImageDifference.getRowDimension(); i++) {
			for(int j = 0; j < Image1AndMeanImageDifference.getColumnDimension()*2; j++)
				if(j < Image1AndMeanImageDifference.getColumnDimension())
					A.setEntry(i, j, Image1AndMeanImageDifference.getEntry(i,j));
				else
					A.setEntry(i, j, Image2AndMeanImageDifference.getEntry(i, j-Image2AndMeanImageDifference.getColumnDimension()));
		}
		System.out.println("A = "+A);
		
		//A transpose times A. ATA.
		RealMatrix ATransposeTimesA = A.transpose().multiply(A);
		//System.out.println("ATransposeTimesA = "+A.transpose()+"*"+A+" = "+ATransposeTimesA);
		System.out.println("ATransposeTimesA = "+ATransposeTimesA);

		//4. Work around: Find the eigenvalues and eigenvectors of A transpose times A.
		EigenDecomposition eigenDecompositionATransposeTimesA = new EigenDecomposition(ATransposeTimesA);
		RealMatrix eigenValueMatrix = eigenDecompositionATransposeTimesA.getD(); //Gets the block diagonal matrix D of the decomposition. D is a block diagonal matrix. Real eigenvalues are on the diagonal while complex values are on 2x2 blocks { {real +imaginary}, {-imaginary, real} }.
		RealMatrix eigenVectorMatrix = eigenDecompositionATransposeTimesA.getV(); //Gets the matrix V of the decomposition. 
		System.out.println("Eigenvalues:");
		System.out.println(eigenValueMatrix);
		System.out.println("Eigenvectors:");
		System.out.println(eigenVectorMatrix);
		
		RealVector[] eigenVector = new RealVector[eigenVectorMatrix.getColumnDimension()];
		for(int i = 0; i < eigenVectorMatrix.getColumnDimension(); i++) {
			eigenVector[i] = eigenDecompositionATransposeTimesA.getEigenvector(i);
		}		
	
		//4.1. Get eigenvalues and eigenvectors of M, not the work around ATA. u = Av.
		RealVector u[] = new RealVector[eigenVector.length];
		for(int i = 0; i < eigenVector.length; i++) {
			u[i] = A.operate(eigenVector[i]); //Eigenvector of M. w = Av.
			System.out.println("Eigenvector of M: "+u[i]);
		}
			
		//5. Compute normalized eigenvector.		 
		RealVector[] normalizedU = new RealVector[u.length];
		for(int i = 0; i < u.length; i++) {			
			/*
			 * double[] magnitude = new double[eigenVector[i].getDimension()]; 
			 * for(int j = 0; j < eigenVector[i].getDimension(); j++) { 
			 * 		magnitude[i] += eigenVector[i].getEntry(j)*eigenVector[i].getEntry(j); 
			 * }
			 * magnitude[i] = Math.sqrt(magnitude[i]); 
			 * normalizedU[i] = eigenVector[i].mapMultiply(1/magnitude[i]); //Normalize: Divide the eigenvector by its own magnitude.
			 */			
			normalizedU[i] = u[i].unitVector();
			System.out.println("Normalized eigenvector: "+normalizedU[i]);
		}
		
		//6. New representation: Representing a face onto the basis of the normalized eigenvector. 
		//I = w*u where w = normalizeUTranspose*ImageAndMeanImageDifference
		
		//normalizedUTranspose
		RealMatrix[] normalizedUTranspose = new RealMatrix[normalizedU.length];		
		for(int i = 0; i < normalizedU.length; i++) {
			normalizedUTranspose[i] = new Array2DRowRealMatrix(normalizedU[i].toArray());
			System.out.println("Normalized eigenvector transpose: "+normalizedUTranspose[i]);
		}
		
//		//A scalar that acts as a weight w_j to the eigenvector in the linear combination expression for representing the original image.
//		RealMatrix I1 = w[0].multiply(Image1AndMeanImageDifference); 
//		RealMatrix I2 = w[1].multiply(Image2AndMeanImageDifference); 
//		
//		for(int i = 0; i < w.length; i++) {
//			w[i] = normalizedUTranspose[i].multiply();
//		}
//		
//		//A scalar that acts as a weight w_j to the eigenvector in the linear combination expression for representing the original image.
//		RealMatrix I1 = w[0].multiply(Image1AndMeanImageDifference); 
//		RealMatrix I2 = w[1].multiply(Image2AndMeanImageDifference); 
//		
//		System.out.println("I1 = "+I1);
//		System.out.println("I2 = "+I2);
//		
//		//7. Feature vector = set of weight = {w1, w2,..., wk} 
//		/*Feature vector = {wTranspose[i]*InputImageAndMeanImageDifference[i]} = {w1, w2,..., wk}.*/
//		
//		//For every image vector, compare vectors.
//		//Compare against a third image.
//		double[][] Image3 = {{9}, {10}, {11}, {12}}; //Centered, normalized, linearized.
//		RealMatrix Image3Matrix = new Array2DRowRealMatrix(Image3);
//		RealMatrix InputImageAndMeanImageDifference3 = Image3Matrix.subtract(meanImage);
//		
//		//Compare wT*InputImageAndMeanImageDifference3 against wT*InputImageAndMeanImageDifference1 and wT*InputImageAndMeanImageDifference2.
//		//Decision function.
	}
}