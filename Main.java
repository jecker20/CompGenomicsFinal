import java.io.FileWriter;
import java.io.IOException;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Random;
import java.util.Hashtable;
import java.util.List;
public class Main {
  public static void main(String[] args) throws IOException {

    //instantiate constant variables

    //bounds of read length
    final int LOWERBOUNDLENGTH = 40;
    final int UPPERBOUNDLENGTH = 50;

    //output file names
    String outputProteinFileName = "ProteinOutputShortTest1.txt";
    String outputReadFileName = "ReadOutputShortTest1.txt";

    //total number of reads
    final int NUMOFREADS = 10;

    //how much of the protein has no read associated
    final double PERCENTOFJUNK = 0;

    //how many of the possible proteins do we store
    final double PERCENTOFPROTEINTOKEEP = 0.2;

    //store reads and amino acids in lists
    List<String[]> reads = new ArrayList<>();
    List<String[]> aminoAcids = new ArrayList<>();

    //file writers
    FileWriter outputRead = new FileWriter(outputReadFileName);
    FileWriter outputProtein = new FileWriter(outputProteinFileName);

    //random vairables
    Random random = new Random();
    String read = "";

    //for loop to create all reads
    for (int i = 0; i < NUMOFREADS; i++) {
      //randomly decides length of read
      int length = random.nextInt((UPPERBOUNDLENGTH + 1) - LOWERBOUNDLENGTH) + LOWERBOUNDLENGTH;

      //makes sure read length of multiple of 3
      length = length - length % 3;

      //creates read
      read = genNewRead(length);

      //stores read and name in array
      String[] newRead = new String[2];
      String readName = "read" + i;
      newRead[0] = readName;
      newRead[1] = read;

      //adds read arrays to a list and prints them to output
      reads.add(newRead);
      outputRead.write(readName + "\n");
      outputRead.write(read + "\n");
    }

    //for loop to iterate over all reading frames of all generated reads
    for (int i = 0; i < reads.size(); i++) {
      //gets read data from array stored in list reads
      String[] readArr = reads.get(i);
      String readName = readArr[0];
      read = readArr[1];
      int length = read.length();

      //creates all reading frames
      String forwardFrame1 = read;
      String forwardFrame2 = forwardFrame1.substring(1, length - 2);
      String forwardFrame3 = forwardFrame1.substring(2, length - 1);
      String backwardFrame1 = reverseRead(new StringBuilder(forwardFrame1).reverse().toString());
      String backwardFrame2 = reverseRead(new StringBuilder(forwardFrame2).reverse().toString());
      String backwardFrame3 = reverseRead(new StringBuilder(forwardFrame3).reverse().toString());

      //converts frames to RNA
      String forwardFrame1RNA = makeRNA(forwardFrame1);
      String forwardFrame2RNA = makeRNA(forwardFrame2);
      String forwardFrame3RNA = makeRNA(forwardFrame3);
      String backwardFrame1RNA = makeRNA(backwardFrame1);
      String backwardFrame2RNA = makeRNA(backwardFrame2);
      String backwardFrame3RNA = makeRNA(backwardFrame3);

      //converts RNA to Amino Acids
      String forwardFrame1AminoAcid = makeAminoAcid(forwardFrame1RNA);
      String forwardFrame2AminoAcid = makeAminoAcid(forwardFrame2RNA);
      String forwardFrame3AminoAcid = makeAminoAcid(forwardFrame3RNA);
      String backwardFrame1AminoAcid = makeAminoAcid(backwardFrame1RNA);
      String backwardFrame2AminoAcid = makeAminoAcid(backwardFrame2RNA);
      String backwardFrame3AminoAcid = makeAminoAcid(backwardFrame3RNA);


      //stores all amino acids in same method as reads
      String[] ff1Arr = new String[2];
      ff1Arr[0] = "read"+i+"-1";
      ff1Arr[1] = forwardFrame1AminoAcid;
      aminoAcids.add(ff1Arr);

      String[] ff2Arr = new String[2];
      ff2Arr[0] = "read"+i+"-2";
      ff2Arr[1] = forwardFrame2AminoAcid;
      aminoAcids.add(ff2Arr);

      String[] ff3Arr = new String[2];
      ff3Arr[0] = "read"+i+"-3";
      ff3Arr[1] = forwardFrame3AminoAcid;
      aminoAcids.add(ff3Arr);

      String[] bf1Arr = new String[2];
      bf1Arr[0] = "read"+i+"-4";
      bf1Arr[1] = backwardFrame1AminoAcid;
      aminoAcids.add(bf1Arr);

      String[] bf2Arr = new String[2];
      bf2Arr[0] = "read"+i+"-5";
      bf2Arr[1] = backwardFrame2AminoAcid;
      aminoAcids.add(bf2Arr);

      String[] bf3Arr = new String[2];
      bf3Arr[0] = "read"+i+"-6";
      bf3Arr[1] = backwardFrame3AminoAcid;
      aminoAcids.add(bf3Arr);

    }

    //for loop to print out amino acids
    for (int i = 0; i<aminoAcids.size(); i++) {
      double rand = random.nextDouble();

      //determines if protein is retained based on set paramaters
      if (rand <= PERCENTOFPROTEINTOKEEP) {
        String[] aminoAcidArr = aminoAcids.get(i);
        String name = aminoAcidArr[0];
        String sequence = aminoAcidArr[1];
        outputProtein.write(">" + name + "\n");
        outputProtein.write(sequence + "\n");
      }
    }

    //for loop to generate junk proteins based on set parameters
    double numJunk = aminoAcids.size() * PERCENTOFJUNK * PERCENTOFPROTEINTOKEEP;
    for (int i = 0; i<numJunk; i++) {
      int randLength = random.nextInt((UPPERBOUNDLENGTH + 1) - LOWERBOUNDLENGTH) + LOWERBOUNDLENGTH;
      String aminoAcid = genAminoAcid(randLength);
      outputProtein.write("readJunk" +i+"\n");
      outputProtein.write(aminoAcid+"" +"\n");
    }
    outputRead.close();
    outputProtein.close();
  }

  //method to provide the complement of a strand (read input should already be reversed)
  public static String reverseRead(String read) {
    String newRead = "";
    for (int i = 0; i < read.length(); i++) {
      char nuc = read.charAt(i);
      if (nuc == 'A') {
        newRead += 'T';
      } else if (nuc == 'C') {
        newRead += 'G';
      } else if (nuc == 'T') {
        newRead += 'A';
      } else if (nuc == 'G') {
        newRead += 'C';
      }
    }
    return newRead;
  }


  //method to convert DNA to RNA
  public static String makeRNA(String read) {
    String newRead = "";
    for (int i = 0; i < read.length(); i++) {
      char nuc = read.charAt(i);
      if (nuc == 'A') {
        newRead += 'U';
      } else if (nuc == 'C') {
        newRead += 'G';
      } else if (nuc == 'T') {
        newRead += 'A';
      } else if (nuc == 'G') {
        newRead += 'C';
      }
    }
    return newRead;
  }

  //method to make RNA into protein
  public static String makeAminoAcid(String read) {
    String aminoAcid = "";
    for (int i = 0; i < read.length(); i += 3) {
      String codon = read.substring(i, i + 3);
      switch (codon) {
        case "GCU":
        case "GCC":
        case "GCA":
        case "GCG":
          aminoAcid += 'A';
          break;
        case "CGU":
        case "CGC":
        case "CGA":
        case "CGG":
        case "AGA":
        case "AGG":
          aminoAcid += 'R';
          break;
        case "AAU":
        case "AAC":
          aminoAcid += 'N';
          break;
        case "GAU":
        case "GAC":
          aminoAcid += 'D';
          break;
        case "UGU":
        case "UGC":
          aminoAcid += 'C';
          break;
        case "CAA":
        case "CAG":
          aminoAcid += 'Q';
          break;
        case "GAA":
        case "GAG":
          aminoAcid += 'E';
          break;
        case "GGU":
        case "GGC":
        case "GGA":
        case "GGG":
          aminoAcid += 'G';
          break;
        case "CAU":
        case "CAC":
          aminoAcid += 'H';
          break;
        case "AUU":
        case "AUC":
        case "AUA":
          aminoAcid += 'I';
          break;
        case "CUU":
        case "CUC":
        case "CUA":
        case "CUG":
        case "UUA":
        case "UUG":
          aminoAcid += 'L';
          break;
        case "AAA":
        case "AAG":
          aminoAcid += 'K';
          break;
        case "AUG":
          aminoAcid += 'M';
          break;
        case "UUU":
        case "UUC":
          aminoAcid += 'F';
          break;
        case "CCU":
        case "CCC":
        case "CCA":
        case "CCG":
          aminoAcid += 'P';
          break;
        case "UCU":
        case "UCC":
        case "UCA":
        case "UCG":
        case "AGU":
        case "AGC":
          aminoAcid += 'S';
          break;
        case "ACU":
        case "ACC":
        case "ACA":
        case "ACG":
          aminoAcid += 'T';
          break;
        case "UGG":
          aminoAcid += 'W';
          break;
        case "UAU":
        case "UAC":
          aminoAcid += 'Y';
          break;
        case "GUU":
        case "GUC":
        case "GUA":
        case "GUG":
          aminoAcid += 'V';
          break;
        default:
          aminoAcid += '-';
      }
    }
    return aminoAcid;
  }

  //method to randomly generate a new read
  public static String genNewRead ( int length){
    String read = "";
    Random random = new Random();
    for (int i = 0; i < length; i++) {
      double val = random.nextDouble();
      if (val >= 0 && val < 0.25) {
        read += 'A';
      } else if (val >= 0.25 && val < 0.5) {
        read += 'T';
      } else if (val >= 0.5 && val < 0.75) {
        read += 'G';
      } else {
        read += 'C';
      }
    }
    return read;
  }

  //method to randomly generate amino acid chain using known frequency values
  public static String genAminoAcid(int length) {
    Random random = new Random();
    String aminoAcid = "";
    for (int i = 0; i < length / 3; i++) {
      int freq = random.nextInt(100);
      if (freq >= 0 && freq < 7.4) {
        aminoAcid += 'A';
      } else if (freq >= 7.4 && freq < 11.6) {
        aminoAcid += 'R';
      } else if (freq >= 11.6 && freq < 16) {
        aminoAcid += 'N';
      } else if (freq >= 16 && freq < 21.9) {
        aminoAcid += 'D';
      } else if (freq >= 21.9 && freq < 25.2) {
        aminoAcid += 'C';
      } else if (freq >= 25.2 && freq < 31) {
        aminoAcid += 'E';
      } else if (freq >= 31 && freq < 34.7) {
        aminoAcid += 'Q';
      } else if (freq >= 34.7 && freq < 42.1) {
        aminoAcid += 'G';
      } else if (freq >= 42.1 && freq < 45) {
        aminoAcid += 'H';
      } else if (freq >= 45 && freq < 48.8) {
        aminoAcid += 'I';
      } else if (freq >= 48.8 && freq < 56.4) {
        aminoAcid += 'L';
      } else if (freq >= 56.4 && freq < 63.6) {
        aminoAcid += 'K';
      } else if (freq >= 63.6 && freq < 65.4) {
        aminoAcid += 'M';
      } else if (freq >= 65.4 && freq < 69.4) {
        aminoAcid += 'F';
      } else if (freq >= 69.4 && freq < 74.4) {
        aminoAcid += 'P';
      } else if (freq >= 74.4 && freq < 82.5) {
        aminoAcid += 'S';
      } else if (freq >= 82.5 && freq < 88.7) {
        aminoAcid += 'T';
      } else if (freq >= 88.7 && freq < 90) {
        aminoAcid += 'W';
      } else if (freq >= 90 && freq < 93.3) {
        aminoAcid += 'Y';
      } else if (freq >= 93.3 && freq < 100) {
        aminoAcid += 'V';
      }
    }
    return aminoAcid;
  }
}