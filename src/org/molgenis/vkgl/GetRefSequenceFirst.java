package org.molgenis.vkgl;

import java.io.*;
import java.net.URL;
import java.util.Scanner;

/**
 * Created by joeri on 2/8/19.
 */
public class GetRefSequenceFirst {


    private static final String MISSING = ".";
    private static int SLEEP = 100;

    public GetRefSequenceFirst(File vkgl_self, File out) throws Exception
    {
        get(vkgl_self, out);
    }

    public void get(File vkgl_self_f, File out) throws Exception
    {
        PrintWriter pw = new PrintWriter(out);

        Scanner vkgl_self = new Scanner(vkgl_self_f);
        String line = vkgl_self.nextLine(); //skip header

        while(vkgl_self.hasNextLine()) {
            line = vkgl_self.nextLine();
            String[] split = line.split("\t", -1);
            if (split.length != 18) {
                throw new Exception("Split lenght must be 18");
            }

            String chrom = split[0];
            String start = split[1];
            String end = split[2];
            String ref = split[3];

            // get the reference genome at this location, add to data, check later
            if (!ref.equals(MISSING))
            {
                int startI = Integer.parseInt(start);
                int endI = Integer.parseInt(end);
                if(endI == 0) { endI = startI; }
                String realhg19bases = getBase(chrom, startI, endI, SLEEP);
                pw.println(ref + "\t" + realhg19bases);
            } else {
                pw.println(ref + "\t" + MISSING);
            }
            pw.flush();
        }
        pw.flush();
        pw.close();

    }

    public static String getBase(String chr, int start, int stop, int sleep) throws Exception
    {
        URL ucsc = new URL("http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr"+chr+":"+start+","+stop);
        BufferedReader getUrlContent = new BufferedReader(new InputStreamReader(ucsc.openStream()));
        String urlLine;
        while ((urlLine = getUrlContent.readLine()) != null)
        {
            // lines starting with < are part of the XML message
            if(!urlLine.startsWith("<"))
            {
                String base = urlLine.toUpperCase();
                getUrlContent.close();
                Thread.sleep(sleep);
                return base;
            }
        }
        getUrlContent.close();
        return MISSING;
    }

    public static void main(String[] args) throws Exception {

        File vkgl_self = new File(args[0]); // /Users/joeri/Projects/VKGL/VKGL_ACMG2015/release_may/VKGL_consensus_may2018_2018-05-31_11_27_44.cleanup.tsv
        File out = new File(args[1]); // /Users/joeri/Projects/VKGL/VKGL_ACMG2015/release_may/VKGL_consensus_hg19check.tsv
        new GetRefSequenceFirst(vkgl_self, out);

    }
}
