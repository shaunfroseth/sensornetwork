import java.util.Arrays;
import java.util.NoSuchElementException;

public class PriorityQueue {
    private Node[] heap;
    private int heapSize, capacity;

    /** Constructor **/
    public PriorityQueue(int capacity) {
        this.capacity = capacity + 1;
        heap = new Node[this.capacity];
        heapSize = 0;
    }

    public boolean isEmpty() {
        return heapSize == 0;
    }

    public int size(){
        return heapSize;
    }

    public void insert(Node n) {
        Node temp = n;

        if (heapSize == capacity - 1){
            grow();
        }

        heap[++heapSize] = temp;

        int i = heapSize;
        while (i != 1 && temp.energy < heap[i / 2].energy) {
            heap[i] = heap[i / 2];
            i /= 2;
        }
        heap[i] = n;
    }

    public Node remove() {
        // NoSuchElement
        if (isEmpty())
            throw new NoSuchElementException("The queue is empty");

        Node returnNode = heap[1];
        Node temp = heap[heapSize--];

        int parent = 1;
        int child = 2;
        while (child <= heapSize) {
            if (child < heapSize && heap[child].energy > heap[child + 1].energy)
                child++;
            if (temp.energy <= heap[child].energy)
                break;

            heap[parent] = heap[child];
            parent = child;
            child *= 2;
        }
        heap[parent] = temp;
        return returnNode;
    }

    private void grow() {
        Node temp[] = new Node[capacity * 2];
        for (int i = 0; i < capacity; i++) {
            temp[i] = heap[i];
        }
        capacity *= 2;
        heap = Arrays.copyOf(heap, capacity);
    }

    @Override
    public String toString() {
        String output = heap[0].toString() + "\n\n";
        for(int i = 1; i < heapSize; i++) {
            output += heap[i] + " \n\n";
        }
        return output;
    }


}

class Node
{
    int node;
    double energy;

    public Node(int node, double energy)
    {
        this.node = node;
        this.energy = energy;
    }

    public String toString()
    {
        return "Node : "+ node +"\nEnergy : "+ energy;
    }
}
