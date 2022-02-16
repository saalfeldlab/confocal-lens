package org.janelia.saalfeldlab.confocallens;

import java.io.IOException;

import ij.ImageJ;
import ij.ImagePlus;
import loci.formats.FormatException;
import loci.plugins.BF;
import ij.IJ;
import ij.ImagePlus;
import ij.io.Opener;

import java.io.File;
import java.util.Arrays;

public class Automation {
	public static void main(String[] args)
	{
		//read lsm files
		String path = "/Users/spreibi/Documents/Janelia/Projects/Saalfeld lense correction/Scope9_20201119_40X_1024X1024/A1_ZB1_T1_BRN_Scope9_20201107_4X4_AF488_AF594_40X_L01.lsm";

		try
		{
			new ImageJ();

			final ImagePlus[] impStack = BF.openImagePlus(path);

			System.out.println( "Number of images: " + impStack.length + " (this should always be 1), we ignore others" );

			if ( impStack.length > 1 )
				throw new RuntimeException( "More than one image was opened, please check the input carefully." );

			final ImagePlus imp = impStack[ 0 ];

			System.out.println(
					"dimensions: " + imp.getStack().getProcessor( 1 ).getWidth() + "x" + imp.getStack().getProcessor( 1 ).getHeight() +
					", channels: " + imp.getNChannels() +
					", z-slices:" + imp.getNSlices() +
					", timepoints: " + imp.getNFrames() );

			imp.resetDisplayRange();
			imp.show();
			
		} catch (FormatException | IOException e) {
			e.printStackTrace();
			System.exit( 0 );
		}

		//generate mips (Stephan: not sure we need that here)
		
		//normalize local contrast brx 127 bry 127 stds 3.0 (all layers)
		
		//montage all layers. least square, translation.
		
		//Align layers. least square (Elastic)
		
		//Auto resize canvas
		
		//Lens correction (All layers)
		
		//Align layers. least square (Elastic)
		
		//output coordinate transform

	}
}
