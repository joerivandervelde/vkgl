package org.molgenis.vkgl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class Vkgl_Acmg2015_Compare {


    public Vkgl_Acmg2015_Compare(File vkgl_self, File vkgl_acmg, File cgd) throws Exception {
        compare(vkgl_self, vkgl_acmg, cgd);
    }

    public static void main(String[] args) throws Exception {

        File vkgl_self = new File(args[0]); // /Users/joeri/Projects/VKGL/VKGL_ACMG2015/release_may/VKGL_consensus_may2018_2018-05-31_11_27_44.cleanup.refcheck.tsv
        File vkgl_acmg = new File(args[1]); // /Users/joeri/Projects/VKGL/VKGL_ACMG2015/release_may/vkgl_april2018_avinput_intervar.hg19_multianno.txt.intervar
        File cgd = new File(args[2]); // /Users/joeri/github/gavin-plus/src/test/resources/bundle_r1.2/CGD_26jun2018.txt.gz
        new Vkgl_Acmg2015_Compare(vkgl_self, vkgl_acmg, cgd);
    }

    public void compare(File vkgl_self_f, File vkgl_acmg_f, File cgd_f) throws Exception {
        // VKGL input data columns expected:
        // chrom	POS	stop	REF	ALT	gene	cDNA	protein	amc	nki	umcg	lumc	vumc	radboud	umcu	erasmus	consensus_classification	disease hg19_ref_check
        //
        // corresponding to numbering used here:
        // 0	    1	2	    3	4	5	    6	    7	    8	9	10	    11	    12	    13	    14	    15	    16	                        17  18
        //
        Scanner vkgl_self = new Scanner(vkgl_self_f);
        int cons_index = 16;

        // InterVar results data columns expected:
        // #Chr	Start	End	Ref	Alt	Ref.Gene	Func.refGene	ExonicFunc.refGene	Gene.ensGene	avsnp147	AAChange.ensGene	AAChange.refGene	clinvar: Clinvar 	 InterVar: InterVar and Evidence 	Freq_ExAC_ALL	Freq_esp6500siv2_all	Freq_1000g2015aug_all	CADD_raw	CADD_phred	SIFT_score	GERP++_RS	phyloP46way_placental	dbscSNV_ADA_SCORE	dbscSNV_RF_SCORE	Interpro_domain	AAChange.knownGene	rmsk	MetaSVM_score	Freq_ExAC_POPs	OMIM	Phenotype_MIM	OrphaNumber	Orpha	Otherinfo
        Scanner vkgl_acmg = new Scanner(vkgl_acmg_f);

        // CHR_START_END_REF_ALT -> gene
        HashMap<String, String> locToGene = new HashMap<>();

        // CHR_START_END_REF_ALT -> B / LB / VUS / LP / P
        HashMap<String, String> locToIntervar = new HashMap<>();

        // CHR_START_END_REF_ALT -> full line of data
        HashMap<String, String> locToIcData = new HashMap<>();


        String line = vkgl_acmg.nextLine(); //skip header
        System.out.println("skipping header: " + line);

        while (vkgl_acmg.hasNextLine()) {
            line = vkgl_acmg.nextLine();
            String[] split = line.split("\t", -1);
            if (split.length != 34) {
                throw new Exception("Split lenght must be 34");
            }

            String chrom = split[0];
            String start = split[1];
            String end = split[2];
            String ref = split[3];
            String alt = split[4];

            String key = chrom + "_" + start + "_" + end + "_" + ref + "_" + alt;

            if (locToGene.containsKey(key)) {
                //  throw new Exception("Duplicate key: " + key);
            }

            locToGene.put(key, split[5]);

            locToIcData.put(key, line);

            String iv = null;
            if (split[13].contains("InterVar: Likely benign")) {
                iv = "B";
            } else if (split[13].contains("InterVar: Likely pathogenic")) {
                iv = "P";
            } else if (split[13].contains("InterVar: Uncertain significance")) {
                iv = "VUS";
            } else if (split[13].contains("InterVar: Benign")) {
                iv = "B";
            } else if (split[13].contains("InterVar: Pathogenic")) {
                iv = "P";
            } else {
                throw new Exception("unknown classification: " + split[13]);
            }

            locToIntervar.put(key, iv);

        }

        System.out.println("done reading intervar output");


        Map<String, CGDEntry> CGD = LoadCGD.loadCGD(cgd_f);

        line = vkgl_self.nextLine(); //skip header
        System.out.println("skipping header: " + line);

        int count = 0;
        int varMatch = 0;
        int refMatch = 0;
        int geneMatch = 0;

        ResultCount allGenes = new ResultCount();

        HashMap<String, ResultCount> countPerGene = new HashMap<>();
        HashMap<String, ArrayList<String>> conflictsPerGene = new HashMap<>();

        while (vkgl_self.hasNextLine()) {
            line = vkgl_self.nextLine();
            String[] split = line.split("\t", -1);
            if (split.length != 19) {
                throw new Exception("Split lenght must be 18");
            }

            String chrom = split[0];
            String start = split[1];
            String end = split[2];
            String ref = split[3];
            String alt = split[4];
            String gene = split[5];
            String ref_hg19check = split[18];


            String key = chrom + "_" + start + "_" + end + "_" + ref + "_" + alt;

            if (locToGene.containsKey(key)) {

                if (ref.equals(ref_hg19check)) {

                    if (locToGene.get(key).equals(gene)) {

                        if (!countPerGene.containsKey(gene)) {
                            countPerGene.put(gene, new ResultCount());
                        }

                        if (!conflictsPerGene.containsKey(gene)) {
                            conflictsPerGene.put(gene, new ArrayList<>());
                        }

                        String consensus = null;
                        if (split[cons_index].contains("pathogenic")) {
                            consensus = "P";
                        } else if (split[cons_index].contains("benign")) {
                            consensus = "B";
                        } else if (split[cons_index].contains("VUS")) {
                            consensus = "VUS";
                        } else if (split[cons_index].contains("Opposite")) {
                            consensus = "CONFLICTING";
                        } else if (split[cons_index].contains("one lab")) {
                            consensus = "ONELAB";
                        } else if (split[cons_index].contains("No consensus")) {
                            consensus = "NOCONS";
                        } //LB+VUS, P+VUS, etc
                        else {
                            throw new Exception("bad consensus status: " + split[cons_index]);
                        }


                        if (consensus.equals("B") && locToIntervar.get(key).equals("B")) {
                            allGenes.ConsB_IvB++;
                            countPerGene.get(gene).ConsB_IvB++;
                        }

                        if (consensus.equals("P") && locToIntervar.get(key).equals("P")) {
                            allGenes.ConsP_IvP++;
                            countPerGene.get(gene).ConsP_IvP++;
                        }

                        if (consensus.equals("B") && locToIntervar.get(key).equals("P")) {
                            allGenes.ConsB_IvP++;
                            countPerGene.get(gene).ConsB_IvP++;
                            conflictsPerGene.get(gene).add(">Conflict: VKGL (consensus): benign, InterVar: pathogenic. Comparison:\n\tVKGL: " + line + "\n\t" + "InterVar: " + locToIcData.get(key));
                        }

                        if (consensus.equals("P") && locToIntervar.get(key).equals("B")) {
                            allGenes.ConsP_IvB++;
                            countPerGene.get(gene).ConsP_IvB++;
                            conflictsPerGene.get(gene).add(">Conflict: VKGL (consensus): pathogenic, InterVar: benign. Comparison:\n\t" + "VKGL: " + line + "\n\t" + "InterVar: " + locToIcData.get(key));
                        }

                        if (consensus.equals("B") && locToIntervar.get(key).equals("VUS")) {
                            allGenes.ConsB_IvVUS++;
                            countPerGene.get(gene).ConsB_IvVUS++;
                        }

                        if (consensus.equals("P") && locToIntervar.get(key).equals("VUS")) {
                            allGenes.ConsP_IvVUS++;
                            countPerGene.get(gene).ConsP_IvVUS++;
                        }

                        if (consensus.equals("VUS") && locToIntervar.get(key).equals("B")) {
                            allGenes.ConsVUS_IvB++;
                            countPerGene.get(gene).ConsVUS_IvB++;
                        }

                        if (consensus.equals("VUS") && locToIntervar.get(key).equals("P")) {
                            allGenes.ConsVUS_IvP++;
                            countPerGene.get(gene).ConsVUS_IvP++;
                        }

                        if (consensus.equals("VUS") && locToIntervar.get(key).equals("VUS")) {
                            allGenes.ConsVUS_IvVUS++;
                            countPerGene.get(gene).ConsVUS_IvVUS++;
                        }

                        if (consensus.equals("NOCONS")) {
                            allGenes.VkglNoCons++;
                            countPerGene.get(gene).VkglNoCons++;
                        }

                        if (consensus.equals("CONFLICTING")) {
                            allGenes.VkglConflict++;
                            countPerGene.get(gene).VkglConflict++;
                        }


                        String oneLabClsf = null;
                        if (consensus.equals("ONELAB")) {
                            String concat = split[8] + split[9] + split[10] + split[11] + split[12] + split[13] + split[14] + split[15];
                            if (concat.contains("Pathogenic") || concat.contains("Likely pathogenic")) {
                                oneLabClsf = "P";
                            } else if (concat.contains("Benign") || concat.contains("Likely benign")) {
                                oneLabClsf = "B";
                            } else if (concat.contains("VUS")) {
                                oneLabClsf = "VUS";
                            } else {
                                throw new Exception("bad one lab classification: " + concat);
                            }

                            if (oneLabClsf.equals("B") && locToIntervar.get(key).equals("B")) {
                                allGenes.OneLabB_IvB++;
                                countPerGene.get(gene).OneLabB_IvB++;
                            }

                            if (oneLabClsf.equals("P") && locToIntervar.get(key).equals("P")) {
                                allGenes.OneLabP_IvP++;
                                countPerGene.get(gene).OneLabP_IvP++;
                            }


                            if (oneLabClsf.equals("B") && locToIntervar.get(key).equals("P")) {
                                allGenes.OneLabB_IvP++;
                                countPerGene.get(gene).OneLabB_IvP++;
                                conflictsPerGene.get(gene).add(">Conflict: VKGL (one lab): benign, InterVar: pathogenic. Comparison:\n\t" + "VKGL: " + line + "\n\t" + "InterVar: " + locToIcData.get(key));
                            }

                            if (oneLabClsf.equals("P") && locToIntervar.get(key).equals("B")) {
                                allGenes.OneLabP_IvB++;
                                countPerGene.get(gene).OneLabP_IvB++;
                                conflictsPerGene.get(gene).add(">Conflict: VKGL (one lab): pathogenic, InterVar: benign. Comparison:\n\t" + "VKGL: " + line + "\n\t" + "InterVar: " + locToIcData.get(key));
                            }

                            if (oneLabClsf.equals("B") && locToIntervar.get(key).equals("VUS")) {
                                allGenes.OneLabB_IvVUS++;
                                countPerGene.get(gene).OneLabB_IvVUS++;
                            }

                            if (oneLabClsf.equals("P") && locToIntervar.get(key).equals("VUS")) {
                                allGenes.OneLabP_IvVUS++;
                                countPerGene.get(gene).OneLabP_IvVUS++;
                            }

                            if (oneLabClsf.equals("VUS") && locToIntervar.get(key).equals("B")) {
                                allGenes.OneLabVUS_IvB++;
                                countPerGene.get(gene).OneLabVUS_IvB++;
                            }

                            if (oneLabClsf.equals("VUS") && locToIntervar.get(key).equals("P")) {
                                allGenes.OneLabVUS_IvP++;
                                countPerGene.get(gene).OneLabVUS_IvP++;
                            }

                            if (oneLabClsf.equals("VUS") && locToIntervar.get(key).equals("VUS")) {
                                allGenes.OneLabVUS_IvVUS++;
                                countPerGene.get(gene).OneLabVUS_IvVUS++;
                            }


                        }
                        geneMatch++;
                    }
                    refMatch++;
                }
                varMatch++;
            }
            count++;

        }

        System.out.println("Gene" + "\t" + ResultCount.header);
        for (String gene : countPerGene.keySet()) {
            System.out.println(gene + "\t" + countPerGene.get(gene).toString());
        }

        System.out.println("TOTAL ALL GENES" + "\t" + allGenes.toString());
        System.out.println("variants seen: " + count + ", variant matched " + varMatch + ", ref matched " + refMatch + ", gene match " + geneMatch);

        // System.exit(1);

        for (String gene : conflictsPerGene.keySet()) {
            if (conflictsPerGene.get(gene).size() == 0) {
                continue;
            }
            System.out.println("\n" + gene + ", LP/P-LB/B conflicts: " + conflictsPerGene.get(gene).size() + " out of " + countPerGene.get(gene).VkglTotal() + "." + (CGD.containsKey(gene) ? " Linked to: " + CGD.get(gene).getCondition() + " (" + CGD.get(gene).getInheritance() + ")" : ""));
            for (String s : conflictsPerGene.get(gene)) {
                System.out.println(s);
            }
        }


        for (String gene : conflictsPerGene.keySet()) {
            if (conflictsPerGene.get(gene).size() > 0) {
                continue;
            }
            System.out.println("\n" + gene + ", LP/P-LB/B conflicts: " + conflictsPerGene.get(gene).size() + " out of " + countPerGene.get(gene).VkglTotal() + "." + (CGD.containsKey(gene) ? " Linked to: " + CGD.get(gene).getCondition() + " (" + CGD.get(gene).getInheritance() + ")" : ""));
        }

    }

}

class ResultCount {
    public int ConsB_IvB = 0;
    public int ConsP_IvP = 0;
    public int ConsB_IvP = 0;
    public int ConsP_IvB = 0;
    public int ConsB_IvVUS = 0;
    public int ConsP_IvVUS = 0;
    public int ConsVUS_IvB = 0;
    public int ConsVUS_IvP = 0;
    public int ConsVUS_IvVUS = 0;
    public int VkglNoCons = 0;
    public int VkglConflict = 0;

    public int OneLabB_IvB = 0;
    public int OneLabP_IvP = 0;
    public int OneLabB_IvP = 0;
    public int OneLabP_IvB = 0;
    public int OneLabB_IvVUS = 0;
    public int OneLabP_IvVUS = 0;
    public int OneLabVUS_IvB = 0;
    public int OneLabVUS_IvP = 0;
    public int OneLabVUS_IvVUS = 0;

    public static String header = "ConsB_IvB" + "\t" + "ConsP_IvP" + "\t" + "ConsB_IvP" + "\t" + "ConsP_IvB" + "\t" + "ConsB_IvVUS" + "\t" + "ConsP_IvVUS" + "\t" + "ConsVUS_IvB" + "\t" + "ConsVUS_IvP" + "\t" + "ConsVUS_IvVUS" + "\t" + "VkglNoCons" + "\t" + "VkglConflict" + "\t" + "OneLabB_IvB" + "\t" + "OneLabP_IvP" + "\t" + "OneLabB_IvP" + "\t" + "OneLabP_IvB" + "\t" + "OneLabB_IvVUS" + "\t" + "OneLabP_IvVUS" + "\t" + "OneLabVUS_IvB" + "\t" + "OneLabVUS_IvP" + "\t" + "OneLabVUS_IvVUS" + "\t" + "ConsTotal" + "\t" + "OneLabTotal" + "\t" + "VkglTotal" + "\t" + "ConflictCons" + "\t" + "ConflictOneLab" + "\t" + "ConflictTotal" + "\t" + "PercConflictCons" + "\t" + "PercConflictOneLab" + "\t" + "PercConflictTotal" + "\t" + "NoConsWithCons" + "\t" + "NoConsWithOneLab" + "\t" + "NoConsWithTotal" + "\t" + "PercNoConsWithCons" + "\t" + "PercNoConsWithOneLab" + "\t" + "PercNoConsWithTotal";

    public int VkglTotal() {
        return ConsB_IvB + ConsP_IvP + ConsB_IvP + ConsP_IvB + ConsB_IvVUS + ConsP_IvVUS + ConsVUS_IvB + ConsVUS_IvP + ConsVUS_IvVUS + VkglNoCons + VkglConflict + OneLabB_IvB + OneLabP_IvP + OneLabB_IvP + OneLabP_IvB + OneLabB_IvVUS + OneLabP_IvVUS + OneLabVUS_IvB + OneLabVUS_IvP + OneLabVUS_IvVUS;
    }

    @Override
    public String toString() {
        int ConsTotal = ConsB_IvB + ConsP_IvP + ConsB_IvP + ConsP_IvB + ConsB_IvVUS + ConsP_IvVUS + ConsVUS_IvB + ConsVUS_IvP + ConsVUS_IvVUS + VkglNoCons + VkglConflict;
        int OneLabTotal = OneLabB_IvB + OneLabP_IvP + OneLabB_IvP + OneLabP_IvB + OneLabB_IvVUS + OneLabP_IvVUS + OneLabVUS_IvB + OneLabVUS_IvP + OneLabVUS_IvVUS;
        int VkglTotal = ConsTotal + OneLabTotal;
        int ConflictCons = ConsB_IvP + ConsP_IvB;
        int ConflictOneLab = OneLabB_IvP + OneLabP_IvB;
        int ConflictTotal = ConflictCons + ConflictOneLab;
        int NoConsWithCons = ConsB_IvVUS + ConsP_IvVUS + ConsVUS_IvB + ConsVUS_IvP;
        int NoConsWithOneLab = OneLabB_IvVUS + OneLabP_IvVUS + OneLabVUS_IvB + OneLabVUS_IvP;
        int NoConsWithTotal = NoConsWithCons + NoConsWithOneLab;

        return ConsB_IvB + "\t" + ConsP_IvP + "\t" + ConsB_IvP + "\t" + ConsP_IvB + "\t" + ConsB_IvVUS + "\t" + ConsP_IvVUS + "\t" + ConsVUS_IvB + "\t" + ConsVUS_IvP + "\t" + ConsVUS_IvVUS + "\t" + VkglNoCons + "\t" + VkglConflict + "\t" + OneLabB_IvB + "\t" + OneLabP_IvP + "\t" + OneLabB_IvP + "\t" + OneLabP_IvB + "\t" + OneLabB_IvVUS + "\t" + OneLabP_IvVUS + "\t" + OneLabVUS_IvB + "\t" + OneLabVUS_IvP + "\t" + OneLabVUS_IvVUS + "\t" + ConsTotal + "\t" + OneLabTotal + "\t" + VkglTotal + "\t" + ConflictCons + "\t" + ConflictOneLab + "\t" + ConflictTotal + "\t" + (ConsTotal == 0 ? 0 : (((double) ConflictCons) / (double) ConsTotal) * 100.0) + "\t" + (OneLabTotal == 0 ? 0 : (((double) ConflictOneLab) / (double) OneLabTotal) * 100.0) + "\t" + (VkglTotal == 0 ? 0 : (((double) ConflictTotal) / (double) VkglTotal) * 100.0) + "\t" + NoConsWithCons + "\t" + NoConsWithOneLab + "\t" + NoConsWithTotal + "\t" + (ConsTotal == 0 ? 0 : (((double) NoConsWithCons) / (double) ConsTotal) * 100.0) + "\t" + (OneLabTotal == 0 ? 0 : (((double) NoConsWithOneLab) / (double) OneLabTotal) * 100.0) + "\t" + (VkglTotal == 0 ? 0 : (((double) NoConsWithTotal) / (double) VkglTotal) * 100.0);
    }
}