package org.molgenis.vkgl;

import org.molgenis.vkgl.CGDEntry.generalizedInheritance;

import java.io.*;
import java.util.Map;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

public class LoadCGD {

	public static Map<String, CGDEntry> loadCGD(File cgdFile) throws IOException
	{

		InputStream fileStream = new FileInputStream(cgdFile);
		InputStream gzipStream = new GZIPInputStream(fileStream);
		Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
		BufferedReader buffered = new BufferedReader(decoder);
		Map<String, CGDEntry> cgd =  new TreeMap<>(String.CASE_INSENSITIVE_ORDER);

		buffered.lines().forEach(line -> {

			String[] split = line.split("\t", -1);
			CGDEntry entry = new CGDEntry(split[0], split[1], split[2], split[3], split[4], split[5], split[6], split[7], split[8], split[9], split[10], split[11]);

			generalizedInheritance inherMode = generalizedInheritance.OTHER;
			if(split[4].contains("AD") && split[4].contains("AR"))
			{
				inherMode = generalizedInheritance.DOMINANT_OR_RECESSIVE;
			}
			else if(split[4].contains("AR"))
			{
				inherMode = generalizedInheritance.RECESSIVE;
			}
			else if(split[4].contains("AD"))
			{
				inherMode = generalizedInheritance.DOMINANT;
			}
			else if(split[4].contains("BG"))
			{
				inherMode = generalizedInheritance.BLOODGROUP;
			}
			else if(split[4].contains("XL"))
			{
				inherMode = generalizedInheritance.X_LINKED;
			}

			entry.setGeneralizedInheritance(inherMode);

			cgd.put(split[0], entry);

		});

		buffered.close();

        return cgd;
	}

	public static void main(String[] args) throws IOException {
		Map<String, CGDEntry> cgd = LoadCGD.loadCGD(new File("/your/path/to/CGD_1jun2016.txt.gz"));

		for( String key : cgd.keySet())
		{
			System.out.println(cgd.get(key).getGene() + " - " + cgd.get(key).getManifestationCategories() + " - " + cgd.get(key).getManifestationCategoriesList().toString());
		}

	}

}
