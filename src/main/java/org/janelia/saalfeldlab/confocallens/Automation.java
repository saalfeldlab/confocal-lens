package org.janelia.saalfeldlab.confocallens;

import java.io.IOException;

import ij.ImageJ;
import ij.ImagePlus;
import loci.formats.FormatException;
import loci.plugins.BF;

public class Automation {
	public static void main(String[] args)
	{
		// test opening
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
	}
}
