/**
 * Exception, thrown when interpolating at a point, which is not within domain, determined by a grid.
 */
public class OutOfGridException extends Exception {

    public OutOfGridException() {}

    public OutOfGridException(String message)
    {
        super(message);
    }

}
