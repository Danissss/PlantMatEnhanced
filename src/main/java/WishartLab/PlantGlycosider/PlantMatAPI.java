package WishartLab.PlantGlycosider;


import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

public class PlantMatAPI {

	/**
	 * generate options from user input
	 * 
	 * @return
	 */
	public static Options generateOptions(){
	
		final Option threadNumber = Option.builder("t")
			.required(false)
			.hasArg(true)
			.argName("Utilizing core")
			.longOpt("thread")
			.desc("Number of core you want to use.")
			.build();
	
	
		final Option filename = Option.builder("f")
			.required(true)
			.hasArg(true)
			.argName("File name")
			.longOpt("file")
			.desc("File name")
			.build();


		final Option filetype = Option.builder("ftype")
			.required(false)
			.hasArg(true)
			.argName("File type")
			.longOpt("ftype")
			.desc("File type: sdf or csv")
			.build();

		final Option output = Option.builder("o")
			.required(true)
			.hasArg(true)
			.argName("Output file")
			.longOpt("output")
			.desc("Define output file location.")
			.build();
		
		final Option csvOutputOption = Option.builder("ocsv")
			.required(false)
			.hasArg(true)
			.argName("Csv Output")
			.longOpt("csvoutput")
			.desc("Select this option to return CSV output(s). You must enter an output filename")
			.build();
	
	
		final Option sdfOutputOption = Option.builder("osdf")
			.required(false)
			.hasArg(true)
			.argName("Sdf Output")
			.longOpt("sdfoutput")
			.desc("Select this option to return SDF output(s). You must enter an output filename")
			.build();

		final Option exception_time = Option.builder("et")
			.required(false)
			.hasArg(true)
			.argName("exception time")
			.longOpt("exceptiontime")
			.desc("The exception time that will terminate the long run thread.")
			.build();
		

		final Option helpOption = Option.builder("h")
			.required(false)
			.hasArg(false)
			.argName("help")
			.longOpt("help")
			.desc("Prints the usage.")
			.build();
	
		final Options options = new Options();
		options.addOption(threadNumber);
		options.addOption(filename);
		options.addOption(filetype);
		options.addOption(output);
		options.addOption(csvOutputOption);
		options.addOption(sdfOutputOption);
		options.addOption(exception_time);
		options.addOption(helpOption);

		return options;
	}
}
