import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;
public class CheckForAccuracy {
  public static void main(String[] args) throws IOException {
    FileReader input = new FileReader("");
    Scanner scan = new Scanner(input);

    int countMatch = 0;
    int countTotal = 0;
    //reads the input file
    while (scan.hasNext()) {
      String line = scan.nextLine();
      int indexOfComma = line.indexOf(',');
      String read = line.substring(0, indexOfComma);
      String protein = line.substring(indexOfComma, line.length());
      if (read.equals(protein)){
        countMatch++;
      }
      countTotal++;
    }

    input.close();

    FileReader input2 = new FileReader("accuracytestprotein.txt");
    int totalProtein = 0;
    scan = new Scanner(input2);
    String readName = scan.nextLine().substring(1,6);
    while (readName.equals("read")) {
      totalProtein++;
    }

    boolean allMatch = (countMatch== countTotal);
    boolean noMissed = (totalProtein == countMatch);

    if (allMatch) {
      System.out.println("All of the read-protein matches are correct");
    } else {
      System.out.println("Some of the read-protein matches are NOT correct");
    }

    if (noMissed) {
      System.out.println("All of the proteins were correctly matched");
    } else {
      System.out.println("Some of the proteins were NOT correctly matched");
    }
  }

}
