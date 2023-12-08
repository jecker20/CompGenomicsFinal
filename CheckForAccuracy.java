import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;
public class CheckForAccuracy {
  public static void main(String[] args) throws IOException {
    //files to be comapred
    String matchedFile = "";
    String proteinFile = "";

    //reads in match file
    FileReader input = new FileReader(matchedFile);
    Scanner scan = new Scanner(input);

    //counter for the number of correct matches
    int countMatch = 0;
    //counter for total read-match pairs
    int countTotal = 0;
    //reads the input file
    while (scan.hasNext()) {
      String line = scan.nextLine();
      int indexOfComma = line.indexOf(',');
      //seperates the read,match line into seperate read and match strings
      String read = line.substring(0, indexOfComma);
      String protein = line.substring(indexOfComma +1, line.length());
      //checks if the read and match are equal and increments appropriately
      if (read.equals(protein)){
        countMatch++;
      }
      countTotal++;
    }

    input.close();

    //opens the protein file
    FileReader input2 = new FileReader(proteinFile);
    int totalProtein = 0;
    scan = new Scanner(input2);
    String readName = scan.nextLine().substring(1,5);
    //counts while the read hasnt reached the junk protein
    while (readName.equals("read") && scan.hasNext()) {
      //increments the total number of expected protein matches
      totalProtein++;
      scan.nextLine();
      if (scan.hasNext()) {
        readName = scan.nextLine().substring(1, 5);
      }
    }

    //checks desired equivalences
    boolean allMatch = (countMatch== countTotal);
    boolean noMissed = (totalProtein == countTotal);

    //ouputs information based on status of equivalences
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
