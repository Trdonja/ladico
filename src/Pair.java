/**
 * A pair of two values, suitable as return type of functions, which return two values.
 */
public class Pair<U, V> {

    public U left;
    public V right;

    public Pair(U left, V right) {
        this.left = left;
        this.right = right;
    }
}
