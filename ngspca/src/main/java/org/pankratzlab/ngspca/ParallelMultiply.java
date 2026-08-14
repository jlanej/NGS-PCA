package org.pankratzlab.ngspca;

import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Matrix multiplication, spread across threads.
 * <p>
 * The products against the full data matrix are what a run at cohort scale spends nearly all of its
 * time on, and {@link BlockRealMatrix#multiply} performs them on one thread.
 * <p>
 * Nothing is reimplemented here. That multiply already walks the output one block-column at a time,
 * so handing it one block-column of the right hand matrix at a time divides exactly the work it was
 * going to do anyway - each output entry still sums its terms in the order it would have, which is
 * why the result is identical to the bit rather than merely close. It is also why the slices are
 * cut on {@link BlockRealMatrix#BLOCK_SIZE}: narrower ones would read the whole left hand matrix
 * once per slice instead of once per block-column, trading arithmetic for memory traffic.
 * <p>
 * That bounds the parallelism at one task per output block-column, so a wider decomposition - a
 * larger {@code numPC} plus oversampling - is what makes more threads useful here.
 */
class ParallelMultiply {

  private ParallelMultiply() {

  }

  /**
   * @param a left hand matrix
   * @param b right hand matrix
   * @param pool the slices are multiplied here
   * @return {@code a} times {@code b}, entry for entry what {@link BlockRealMatrix#multiply} gives
   */
  static RealMatrix multiply(BlockRealMatrix a, RealMatrix b, ForkJoinPool pool) {
    int columns = b.getColumnDimension();
    if (pool.getParallelism() < 2 || columns <= BlockRealMatrix.BLOCK_SIZE) {
      return a.multiply(b);
    }
    int rows = b.getRowDimension();
    int slices = (columns + BlockRealMatrix.BLOCK_SIZE - 1) / BlockRealMatrix.BLOCK_SIZE;

    RealMatrix[] parts = new RealMatrix[slices];
    pool.submit(() -> IntStream.range(0, slices).parallel().forEach(slice -> {
      int from = slice * BlockRealMatrix.BLOCK_SIZE;
      int to = Math.min(from + BlockRealMatrix.BLOCK_SIZE, columns) - 1;
      parts[slice] = a.multiply(b.getSubMatrix(0, rows - 1, from, to));
    })).join();

    BlockRealMatrix product = new BlockRealMatrix(a.getRowDimension(), columns);
    for (int slice = 0; slice < slices; slice++) {
      int from = slice * BlockRealMatrix.BLOCK_SIZE;
      for (int column = 0; column < parts[slice].getColumnDimension(); column++) {
        product.setColumn(from + column, parts[slice].getColumn(column));
      }
      parts[slice] = null;
    }
    return product;
  }
}
