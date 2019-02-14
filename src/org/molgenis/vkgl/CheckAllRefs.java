package org.molgenis.vkgl;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URL;
import java.util.Scanner;

/**
 * Created by joeri on 2/8/19.
 */
public class CheckAllRefs {


    private static final String MISSING = ".";
    private static int SLEEP = 100;

    public CheckAllRefs(File vkgl_self, File out) throws Exception
    {
        get(vkgl_self, out);
    }

    public void get(File vkgl_self_f, File out) throws Exception
    {
        PrintWriter pw = new PrintWriter(out);

        Scanner vkgl_self = new Scanner(vkgl_self_f);
        String line = vkgl_self.nextLine(); //skip header

        pw.println("chrom" + "\t" + "start" + "\t" + "end"  + "\t" + "ref" + "\t" + "hg19refcheck");

        while(vkgl_self.hasNextLine()) {
            line = vkgl_self.nextLine();
            String[] split = line.split("\t", -1);
            if (split.length != 19) {
                throw new Exception("Split lenght must be 19");
            }

            String chrom = split[0];
            String start = split[1];
            String end = split[2];
            String ref = split[3];
            String hg19ref = split[18];

            if(!ref.equals(hg19ref))
            {
                pw.println(chrom + "\t" + start + "\t" + end + "\t" + ref + "\t" + hg19ref);
            }

            pw.flush();
        }
        pw.flush();
        pw.close();

    }

    public static void main(String[] args) throws Exception {

        File vkgl_self = new File(args[0]); // /Users/joeri/Projects/VKGL/VKGL_ACMG2015/release_may/VKGL_consensus_may2018_2018-05-31_11_27_44.cleanup.tsv
        File out = new File(args[1]); // /Users/joeri/Projects/VKGL/VKGL_ACMG2015/release_may/VKGL_consensus_hg19check.tsv
        new CheckAllRefs(vkgl_self, out);

    }
}
