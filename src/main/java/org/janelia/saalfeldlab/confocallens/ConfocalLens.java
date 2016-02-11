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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.collections.MapUtils;
import org.apache.commons.io.filefilter.AndFileFilter;
import org.apache.commons.io.filefilter.FileFileFilter;
import org.apache.commons.io.filefilter.RegexFileFilter;
import org.imagearchive.lsm.reader.Reader;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.Opener;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import mpicbg.models.CoordinateTransform;
import mpicbg.trakem2.transform.NonLinearCoordinateTransform;
import mpicbg.trakem2.transform.TransformMesh;
import mpicbg.trakem2.transform.TransformMeshMappingWithMasks;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccessible;
import net.imglib2.converter.Converter;
import net.imglib2.converter.Converters;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.DoubleArray;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgs;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.RealViews;
import net.imglib2.realtransform.Translation3D;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.ARGBType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.view.Views;

/**
 *
 *
 * @author Stephan Saalfeld <saalfelds@janelia.hhmi.org>
 */
public class ConfocalLens {

	final static public String[] ls(final String dirStr, final String regex) {
		final File dir = new File(dirStr);
		if (dir.exists() && dir.isDirectory()) {
			final String[] fileNames = dir.list(
					new AndFileFilter(
							FileFileFilter.FILE,
							new RegexFileFilter(regex)));
			Arrays.sort(fileNames);
			return fileNames;

		}
		else return new String[0];
	}

	final static public HashMap<String, Double> parseZOffset(final String tileConfigurationStr) throws IOException {

		final HashMap<String, Double> offsets = new HashMap<String, Double>();

		try (
				final BufferedReader reader =
				new BufferedReader(
						new FileReader(
								tileConfigurationStr))) {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				if (line.matches("^[^;]*;.*\\(.*\\)$")) {
					final String fileStr = line.replaceAll("^([^;]*);.*$", "$1");
					final String zStr = line.replaceAll("^[^;]*;.*\\([^,]*,\\s*[^,]*,\\s*([^\\s]*)\\)$", "$1");
					final Double z = Double.parseDouble(zStr);
					offsets.put(fileStr, z);
				}
			}
		}

		return offsets;
	}

	final static public ImagePlus openImagePlus(
			final String dirStr,
			final String fileStr) {
		final ImagePlus imp;
		if (fileStr.endsWith(".tif"))
			imp = new Opener().openImage(dirStr + fileStr);
		else if (fileStr.endsWith(".lsm"))
			imp = new Reader().open(dirStr + fileStr);
		else
			imp = null;

		return imp;
	}

	final static public <T extends NumericType<T> & NativeType<T>> ImagePlusImg<T, ?> openStack(
			final String dirStr,
			final String fileStr) {
		System.out.println(dirStr + fileStr);
		final ImagePlus imp = openImagePlus(dirStr, fileStr);
		if (imp != null)
			return ImagePlusImgs.from(imp);
		else
			return null;
	}

	final static public <T extends NumericType<T> & NativeType<T>> RandomAccessibleInterval<T> zShift(
			final RandomAccessibleInterval<T> src, final double zShift) {
		final RandomAccessible<T> extended = Views.extendBorder(src);
		final RealRandomAccessible<T> interpolant = Views.interpolate(extended, new NLinearInterpolatorFactory<T>());
		final Translation3D shift = new Translation3D(0, 0, zShift);
		final RandomAccessible<T> shiftedRaster = RealViews.affine(interpolant, shift);
		return Views.interval(shiftedRaster, src);
	}

	@SuppressWarnings("unchecked")
	final static public <T extends NumericType<T> & NativeType<T>> FloatProcessor createThickSlice(
			final RandomAccessibleInterval<T> src,
			final int z,
			final int radius) {

		final RandomAccessibleInterval<? extends RealType<?>> realSrc;
		final T t = src.randomAccess().get();
		if (RealType.class.isInstance(t)) {
			realSrc = (RandomAccessibleInterval<RealType<?>>)src;
		} else if (ARGBType.class.isInstance(t)) {
			realSrc = Converters.convert(
					(RandomAccessibleInterval<ARGBType>)src,
					new Converter<ARGBType, DoubleType>(){

						@Override
						public void convert(final ARGBType input, final DoubleType output) {
							final int argb = input.get();

							final double r = (argb >> 16) & 0xff;
							final double g = (argb >> 8) & 0xff;
							final double b = argb & 0xff;

							output.setReal(.3 * r + .6 * g + .1 * b);
						}
					},
					new DoubleType());
		} else {
			return null;
		}

		/* accumulate */
		final long firstZ = Math.max(0, Math.min(src.max(2), z - radius));
		final long lastZ = Math.max(0, Math.min(src.max(2), z + radius));
		final ArrayImg<DoubleType, DoubleArray> acc = ArrayImgs.doubles(src.dimension(0), src.dimension(1));
		for (long zi = firstZ; zi <= lastZ; ++zi) {
			final RandomAccessibleInterval<? extends RealType<?>> slice = Views.hyperSlice(realSrc, 2, zi);
			final IterableInterval<? extends RealType<?>> iterableSlice = Views.flatIterable(slice);
			final Cursor<? extends RealType<?>> iterableSliceCursor = iterableSlice.cursor();
			final Cursor<DoubleType> iterableAccCursor = acc.cursor();
			while (iterableAccCursor.hasNext()) {
				iterableSliceCursor.fwd();
				iterableAccCursor.fwd();
				final DoubleType a = iterableAccCursor.get();
				a.set(a.get() + iterableSliceCursor.get().getRealDouble());
			}
		}

		/* normalize */
		final long n = lastZ - firstZ + 1;
		final Cursor<DoubleType> iterableAccCursor = acc.cursor();
		while (iterableAccCursor.hasNext()) {
			iterableAccCursor.fwd();
			final DoubleType a = iterableAccCursor.get();
			a.set(a.get() / n);
		}

		/* copy to float */
		final float[] floats = new float[(int)src.dimension(0) * (int)src.dimension(1)];
		final Cursor<DoubleType> doubleCursor = acc.cursor();
		for (int i = 0; doubleCursor.hasNext(); ++i) {
			floats[i] = doubleCursor.next().getRealFloat();
		}

		return new FloatProcessor((int)src.dimension(0), (int)src.dimension(1), floats);

	}

	final static public ImageStack createTransformedStack(
			final ImageStack srcStack,
			final CoordinateTransform t,
			final int cropWidth,
			final int meshResolution) {
		final TransformMesh mesh = new TransformMesh(t, 128, srcStack.getWidth(), srcStack.getHeight());
		final Rectangle bounds = mesh.getBoundingBox();
		final ImageStack stack = new ImageStack(bounds.width - 2 * cropWidth, bounds.height - 2 * cropWidth);
		final TransformMeshMappingWithMasks<TransformMesh> mapping = new TransformMeshMappingWithMasks<TransformMesh>(mesh);
		for (int i = 0; i < srcStack.getSize(); ++i) {
			final ImageProcessor src = srcStack.getProcessor(i + 1);
			src.setInterpolationMethod(ImageProcessor.BILINEAR);
			final ImageProcessor dst = src.createProcessor(bounds.width, bounds.height);
			mapping.mapInterpolated(src, dst);
			final ImageProcessor crop = dst.createProcessor(bounds.width - 2 * cropWidth, bounds.height - 2 * cropWidth);
			crop.copyBits(dst, -cropWidth, -cropWidth, Blitter.COPY);

			stack.addSlice(crop);
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
			final String dirStr,
			final String fileName,
			final CoordinateTransform t,
			final int cropWidth) {
		final ImagePlus imp = openImagePlus(dirStr, fileName);
		if (imp != null) {
			imp.setStack(createTransformedStack(imp.getStack(), t, cropWidth));
			return imp;
		}
		return null;
	}

	final static public void showTransformedImages(
			final String dirStr,
			final Iterable<String> fileNames,
			final CoordinateTransform t,
			final int cropWidth) {
		for (final String fileName : fileNames) {
			final ImagePlus imp = loadAndTransformImagePlus(
					dirStr,
					fileName,
					t,
					cropWidth);
			if (imp != null) {
				imp.show();
			}
		}
	}

	@SuppressWarnings({"unchecked", "rawtypes"})
	final static public void createThickSliceStacks(final String dirStr, final Map<String, Double> offsets, final int stepSize) {

		for (int z = 50; z < 350; z += stepSize)
		{
			ImageStack stack = null;
			boolean first = true;
			ImagePlus imp = null;
			for (final Entry<String, Double> entry : offsets.entrySet()) {
				final FloatProcessor ip = createThickSlice(
						(RandomAccessibleInterval)zShift(
								(RandomAccessibleInterval)openStack(
										dirStr,
										entry.getKey()),
								entry.getValue().doubleValue()), z, stepSize / 2);
				if (first)
					stack = new ImageStack(ip.getWidth(), ip.getHeight());

				stack.addSlice(entry.getKey(), ip);
				if (first) {
					imp = new ImagePlus("" + z, stack);
					imp.show();
					first = false;
				}
				imp.setStack(stack);
				imp.updateAndDraw();
			}
		}
	}



	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(final String[] args) throws IOException {
		new ImageJ();

//		final String dirStr = "/groups/saalfeld/home/saalfelds/experiments/fly-light-lens/20151214_SS_63X_brain_tiles/data_1/";
//		final String dirStr = "/groups/saalfeld/home/saalfelds/experiments/fly-light-lens/chem tag 4X4 tile mosaic_scope3_60 percent overlap/data_1/";
//		final String dirStr = "/groups/saalfeld/home/saalfelds/experiments/fly-light-lens/production data_scope 3_2 brains_5-tiles_63X/JRC_SS03173-20150731_21_A3~63x/";
//		final String dirStr = "/groups/saalfeld/home/saalfelds/experiments/fly-light-lens/production data_scope 3_2 brains_5-tiles_63X/JRC_SS04421-20151124_32_J2~63x/";
		final String dirStr = "/misc/public/YoshiToStephan/DPX_Standard2013/lsm_files/lens-model-scope3/";

//		final HashMap<String, Double> offsets = parseZOffset(dirStr + "TileConfiguration.registered.txt");
		final HashMap<String, Double> offsets = parseZOffset(dirStr + "TileConfiguration.txt");

		MapUtils.debugPrint(System.out, "offsets", offsets);

		createThickSliceStacks(dirStr, offsets, 50);


		/* warp all the stacks */
		final NonLinearCoordinateTransform t = new NonLinearCoordinateTransform();

		/* scope X, chemical, z=250 */
		t.init("5 21 265.3552814678547 -9.46667068228438 -4.531760742006943 239.37726894197888 22.858328199010067 12.577412433666721 9.975689928417305 16.35325129865469 9.650389592551617 36.16151963840831 -17.244014738564207 -10.634169785166533 -7.298201392062338 -15.37001738770205 -11.616095930236812 -2.651601920760191 -11.547677405817334 -51.422129456488506 4.2284172035821666 14.066659062796418 4.574271409908075 2.0346108475675018 6.327703558886167 -0.1867516216106493 4.66024082852573 3.098814637577937 9.435288465289805 38.46100321049568 -2.2137394043142162 -6.312635013840618 -0.5157633863136688 -1.4507276842512062 -2.7121668017568448 1.5113537249434053 -0.9037087898330007 -1.3848042937070646 -1.857603111460449 -0.5318306705492479 -2.6272063051772587 -14.671304256902848 5.503570021583938 6.5815444429589665 550.3611303987366 658.1659554903936 378587.47791847424 353367.6744855716 496668.64115102845 2.878994708513267E8 2.3964678240518382E8 2.6408401718693027E8 4.0013589260626554E8 2.3241774823386978E11 1.8063404955169467E11 1.7776898337858954E11 2.117299481573544E11 3.360791915126304E11 1.9539835409852516E14 1.4493774500636447E14 1.3324672546785247E14 1.4183003173533994E14 1.7734286817446294E14 2.9067462900420606E14 100.0 275.12560452535257 251.97136023872915 298486.1321366977 230006.6555085441 299005.3010550548 2.9372703781410813E8 2.2138110056354192E8 2.2053311689287165E8 3.1111484269634515E8 2.840251131113193E11 2.090331229438189E11 1.9579593692097876E11 2.1195446019972284E11 3.115882961767197E11 2.7379440754202656E14 1.9784325948288825E14 1.7858717429114544E14 1.8015632128486216E14 2.037798417056104E14 3.074332037689868E14 0.0 1024 1024 ");

		/* scope 3, chemical, z=150 */
//		t.init("5 21 265.6307143127233 -9.130056188368174 -5.042534502175269 285.04252007260436 16.73750463611927 8.687376653986409 10.115997261395686 15.420559476744526 9.023895052786482 22.821906826633317 -11.529566512492051 1.4000481520338184 -7.246043792423219 -13.939336757579941 -11.91361713152407 -3.1591756660950754 -8.641052066680258 -22.680446946664638 -0.7356219693276049 -2.196917800709194 3.375523369800531 0.4539356754922039 7.612048812901785 2.2811499578695624 4.595327352940267 1.8899066385111287 6.801776782400479 11.777216765250174 1.7810157593166003 1.0785643897196024 0.6865956787357633 -0.39895368665555164 -4.120482963135277 0.3081569949647477 -0.48676918230487964 -1.6517350116837806 -1.9910132841425692 0.1719597431666795 -1.781006585768337 -5.039542303265364 3.98118187165964 6.084927267255223 398.1132481342764 608.4999799578544 233377.64113242135 243931.7229618072 457591.35727234953 1.6407630730792004E8 1.4431750089423496E8 1.839906145493179E8 3.738953280525476E8 1.2789121253936269E11 1.0215596658039322E11 1.0935180677452296E11 1.505606297400481E11 3.1923073603342224E11 1.0611278880000802E14 7.993067470049967E13 7.76977922690329E13 8.967036665094672E13 1.2863916133117648E14 2.801306764577319E14 100.0 273.66777413450376 295.51895782047615 270991.6074233788 223284.92152850868 331445.9034230781 2.5433458201219186E8 2.0031821110479888E8 2.069616917160562E8 3.32811174009026E8 2.386302151222804E11 1.8273372946067606E11 1.7584555817430225E11 1.9374160646313586E11 3.265268357885264E11 2.253879526293362E14 1.6916041762775225E14 1.5750945002548253E14 1.599289099141636E14 1.829640463375873E14 3.186985651267883E14 0.0 1024 1024 ");
//		showTransformedImages(dirStr, offsets.keySet(), t, 16);
	}

}
