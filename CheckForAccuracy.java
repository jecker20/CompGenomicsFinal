import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;
public class CheckForAccuracy {
  public static void main(String[] args) throws IOException {
    FileReader input = new FileReader("out (3).txt");
    Scanner scan = new Scanner(input);

    int countMatch = 0;
    int countTotal = 0;
    //reads the input file
    while (scan.hasNext()) {
      String line = scan.nextLine();
      int indexOfComma = line.indexOf(',');
      String read = line.substring(0, indexOfComma);
      String protein = line.substring(indexOfComma +1, line.length());
      if (read.equals(protein)){
        countMatch++;
      }
      countTotal++;
    }

    input.close();

    FileReader input2 = new FileReader("smallaccuracytestprotein2.txt");
    int totalProtein = 0;
    scan = new Scanner(input2);
    String readName = scan.nextLine().substring(1,5);
    while (readName.equals("read") && scan.hasNext()) {
      totalProtein++;
      scan.nextLine();
      if (scan.hasNext()) {
        readName = scan.nextLine().substring(1, 5);
      }
    }

    boolean allMatch = (countMatch== countTotal);
    boolean noMissed = (totalProtein == countTotal);

    if (allMatch) {
      System.out.println("Read-Protein Matches: NO ERRORS");
    } else {
      System.out.println("Read-Protein Matches: ERRORS");
    }

    if (noMissed) {
      System.out.println("Total Protein: NO ERRORS");
    } else {
      System.out.println("Total Protein: ERRORS");
    }
  }

}
