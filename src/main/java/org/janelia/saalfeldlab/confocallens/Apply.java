/**
 * License: GPL
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 2
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.janelia.saalfeldlab.confocallens;

import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;

import org.imagearchive.lsm.reader.Reader;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.Opener;
import ij.process.Blitter;
import ij.process.ImageProcessor;
import mpicbg.models.CoordinateTransform;
import mpicbg.trakem2.transform.NonLinearCoordinateTransform;
import mpicbg.trakem2.transform.TransformMesh;
import mpicbg.trakem2.transform.TransformMeshMappingWithMasks;

/**
 * Apply a 2D transformation to all slices of a stack that is either lsm or tif.
 * Save as tif.
 *
 * @author Stephan Saalfeld <saalfelds@janelia.hhmi.org>
 */
public class Apply {

	final static public ImagePlus openImagePlus(
			final String path) {
		final ImagePlus imp;
		if (path.endsWith(".tif"))
			imp = new Opener().openImage(path);
		else if (path.endsWith(".lsm"))
			imp = new Reader().open(path);
		else
			imp = null;

		return imp;
	}

	final static public ImageStack createTransformedStack(
			final ImageStack srcStack,
			final CoordinateTransform t,
			final int cropWidth,
			final int meshResolution) {
		final TransformMesh mesh = new TransformMesh(t, meshResolution, srcStack.getWidth(), srcStack.getHeight());
		final Rectangle bounds = mesh.getBoundingBox();
		final ImageStack stack = new ImageStack(bounds.width - 2 * cropWidth, bounds.height - 2 * cropWidth);
		final TransformMeshMappingWithMasks<TransformMesh> mapping = new TransformMeshMappingWithMasks<TransformMesh>(mesh);
		for (int i = 0; i < srcStack.getSize(); ++i) {
			final ImageProcessor src = srcStack.getProcessor(i + 1);
			src.setInterpolationMethod(ImageProcessor.BILINEAR);
			final ImageProcessor dst = src.createProcessor(bounds.width, bounds.height);
			mapping.mapInterpolated(src, dst);
			final ImageProcessor cropped = dst.createProcessor(bounds.width - 2 * cropWidth, bounds.height - 2 * cropWidth);
			cropped.copyBits(dst, -cropWidth, -cropWidth, Blitter.COPY);

			stack.addSlice(cropped);
		}
		return stack;
	}

	final static public ImageStack createTransformedStack(
			final ImageStack srcStack,
			final CoordinateTransform t,
			final int cropWidth) {
		return createTransformedStack(srcStack, t, cropWidth, 128);
	}



	/**
	 * Load an {@link ImagePlus}, transform all its slices with a
	 * {@link CoordinateTransform}, and crop its borders.
	 *
	 * @param dirStr
	 * @param fileName
	 * @param t
	 * @param cropWidth
	 * @return
	 */
	final static public ImagePlus loadAndTransformImagePlus(
			final String path,
			final CoordinateTransform t,
			final int cropWidth) {
		final ImagePlus imp = openImagePlus(path);
		if (imp != null) {
			imp.setStack(createTransformedStack(imp.getStack(), t, cropWidth));
			return imp;
		}
		return null;
	}

	final static public void saveTransformedImages(
			final String dirStr,
			final Iterable<String> fileNames,
			final String outDirStr,
			final CoordinateTransform t,
			final int cropWidth) {
		for (final String fileName : fileNames) {
			final ImagePlus imp = loadAndTransformImagePlus(
					dirStr + fileName,
					t,
					cropWidth);
			if (imp != null) {
				IJ.saveAsTiff(imp, outDirStr + fileName + ".tif");
			}
		}
	}

	static private ImagePlus impInput = null;
	static private String pathOutput = null;
	static private CoordinateTransform transform;
	static private int crop = 0;

	static public boolean setup(final String... args) {
		if (args.length < 4) return false;
		impInput = openImagePlus(args[0]);
		if (impInput == null) return false;
		pathOutput = args[1];
		try {
			final File f = new File(pathOutput).getParentFile();
			if (f == null || !(f.mkdirs() || f.exists())) return false;
			System.out.println(f);
		}
		catch (final Exception e) {
			e.printStackTrace(System.err);;
			return false;
		}
		final NonLinearCoordinateTransform t = new NonLinearCoordinateTransform();
		t.init(args[2]);
		transform = t;
		crop = Integer.parseInt(args[3]);
		return true;
	}


	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(final String[] args) throws IOException {

		if (setup(args)) {
			impInput.setStack(createTransformedStack(impInput.getStack(), transform, crop));
			IJ.saveAsTiff(impInput, pathOutput );
		} else {
			System.err.println("Usage: "
					+ "java ... <input_path> <output_path> \"<lens_model>\" <crop_width>."
					);
		}
	}
}
