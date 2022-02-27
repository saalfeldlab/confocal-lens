package org.janelia.saalfeldlab.confocallens;

import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.FileUtils;
import org.json.JSONObject;
import org.json.JSONStringer;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.HyperStackConverter;
import ij.plugin.ZProjector;
import ij.process.Blitter;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ini.trakem2.ControlWindow;
import ini.trakem2.Project;
import ini.trakem2.display.Display;
import ini.trakem2.display.Displayable;
import ini.trakem2.display.Layer;
import ini.trakem2.display.LayerSet;
import ini.trakem2.display.Patch;
import lenscorrection.DistortionCorrectionTask;
import lenscorrection.DistortionCorrectionTask.CorrectDistortionFromSelectionParam;
import loci.plugins.BF;
import mpicbg.ij.plugin.NormalizeLocalContrast;

import mpicbg.models.*;
//import mpicbg.models.AffineModel2D;
//import mpicbg.models.CoordinateTransformList;
//import mpicbg.models.IdentityModel;
//import mpicbg.models.IllDefinedDataPointsException;
//import mpicbg.models.Model;
//import mpicbg.models.NotEnoughDataPointsException;
//import mpicbg.models.Point;
//import mpicbg.models.PointMatch;
//import mpicbg.models.RigidModel2D;

import mpicbg.trakem2.align.Align;
import mpicbg.trakem2.align.AlignTask;
import mpicbg.trakem2.align.RegularizedAffineLayerAlignment;
import mpicbg.trakem2.transform.CoordinateTransform;

import java.awt.geom.AffineTransform;
import java.awt.image.IndexColorModel;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class Automation {

	static String[] lambdas = new String[]{
		"488nm, pass1",
		"594nm",
		"488nm, pass2",
		"561nm",
		"647nm"};

	static String format = "%s, %s, %s";
	
	static Class<?> invarianceModelClass = IdentityModel.class;
	
	static ImageProcessor visualizeDifference(
			int w,
			int h,
			int pWidth,
			int pHeight,
			CoordinateTransform ct1,
			CoordinateTransform ct2) {
		double sx = (double)pWidth / w;
		double sy = (double)pHeight / h;
		FloatProcessor ip = new FloatProcessor(w, h);
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				double[] l1 = new double[]{x * sx, y * sy};
				double[] l2 = new double[]{x * sx, y * sy};
				ct1.applyInPlace(l1);
				ct2.applyInPlace(l2);
				double dx = l1[0] - l2[0];
				double dy = l1[1] - l2[1];
				double d = Math.sqrt(dx * dx + dy * dy);
				ip.setf(x, y, (float)d);
			}
		}
		//matrix.copyBits( ip, i * 256, j * 256, Blitter.COPY );
		return ip;
	}

	static ImageProcessor visualizeDifferenceVectorDistribution(
			int w,
			int h,
			int pWidth,
			int pHeight,
			CoordinateTransformList<?> ct1,
			CoordinateTransformList<?> ct2,
			double max) {
		double hw = 0.5 * w;
		double hh = 0.5 * h;
		double sx = (double)pWidth / w;
		double sy = (double)pHeight / h;
		FloatProcessor ip = new FloatProcessor(w, h);
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				double[] l1 = new double[]{x * sx, y * sy};
				double[] l2 = new double[]{x * sx, y * sy};
				ct1.applyInPlace(l1);
				ct2.applyInPlace(l2);
				double dx = l1[0] - l2[0];
				double dy = l1[1] - l2[1];
				dx = Math.min(w - 1, Math.max(0.0, (dx / max + 1) * hw));
				dy = Math.min(w - 1, Math.max(0.0, (dy / max + 1) * hh));
				int ix = (int)Math.round(dx);
				int iy = (int)Math.round(dy);
				ip.setf(ix, iy, (float)(1.0 + ip.getPixelValue(ix, iy)));
			}
		}
		//matrix.copyBits( ip, i * 256, j * 256, Blitter.COPY );
		return ip;
	}

	static ImageProcessor visualizeDifferenceVectors(
			int w,
			int h,
			int pWidth,
			int pHeight,
			CoordinateTransformList<?> ct1,
			CoordinateTransformList<?> ct2,
			double max) {
		double sx = (double)pWidth / w;
		double sy = (double)pHeight / h;
		ColorProcessor ip = new ColorProcessor(w, h);
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				double[] l1 = new double[]{x * sx, y * sy};
				double[] l2 = new double[]{x * sx, y * sy};
				ct1.applyInPlace(l1);
				ct2.applyInPlace(l2);
				double dx = (l1[0] - l2[0] ) / max;
				double dy = (l1[1] - l2[1] ) / max;
				double d = Math.sqrt(dx * dx + dy * dy);
				double s = 1.0 / d;
				if (s < 1.0) {
					dx *= s;
					dy *= s;
				}
				
				ip.set(x, y, mpicbg.ij.util.Util.colorVector(dx, dy));
			}
		}
		return ip;
	}

	static Model<?> sampleModel(CoordinateTransform ct, Class<?> modelClass, int width, int height) throws NotEnoughDataPointsException, IllDefinedDataPointsException, InstantiationException, IllegalAccessException {
		Model<?> model = (Model<?>) modelClass.newInstance();
		ArrayList<PointMatch> matches = new ArrayList<PointMatch>();
		double scaleX = ((double)width - 1.0f) / 63.0f;
		double scaleY = ((double)height - 1.0f) / 63.0f;
		for (int y = 0; y < 64; ++y) {
			double ys = scaleY * y;
			for (int x = 0; x < 64; ++x) {
				double xs = scaleX * x;
				Point p = new Point(new double[]{xs, ys});
				p.apply(ct);
				matches.add(new PointMatch(p, p));
			}
		}
		model.fit(matches);
		return model;
	}

	/**
	 * Estimate a transformation model between two CoordinateTransforms.
	 * The model maps ct1 into ct2.
	 * @throws IllDefinedDataPointsException 
	 * @throws NotEnoughDataPointsException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 */
	static Model<?> sampleModel2(
			CoordinateTransformList<?> ct1,
			CoordinateTransformList<?> ct2,
			Class<?> modelClass,
			int width,
			int height) throws NotEnoughDataPointsException, IllDefinedDataPointsException, InstantiationException, IllegalAccessException {
		Model<?> model = (Model<?>) modelClass.newInstance();
		ArrayList<PointMatch> matches = new ArrayList<PointMatch>();
		double scaleX = ((double)width - 1.0f) / 63.0f;
		double scaleY = ((double)height - 1.0f) / 63.0f;
		for (int y = 0; y < 64; ++y) {
			double ys = scaleY * y;
			for (int x = 0; x < 64; ++x) {
				double xs = scaleX * x;
				Point p = new Point(new double[]{xs, ys});
				Point q = new Point(new double[]{xs, ys});
				ct1.applyInPlace(p.getL());
				q.apply(ct2);
				matches.add(new PointMatch(p, q));
			}
		}
		model.fit(matches);
		return model;
	}


	static CoordinateTransform createTransform(
			String className,
			String dataString) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		CoordinateTransform ct = (CoordinateTransform)Class.forName(className).newInstance();
		ct.init(dataString);
		return ct;
	}

	static CoordinateTransformList<mpicbg.models.CoordinateTransform> createTransformList(int j, String[][] transforms) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		CoordinateTransformList<mpicbg.models.CoordinateTransform> ctl = new CoordinateTransformList<mpicbg.models.CoordinateTransform>();
		for (int k = 2; k < transforms[j].length; k +=2)
			ctl.add(createTransform(transforms[j][k-1], transforms[j][k]));
		return ctl;
	}

	static ImagePlus showDifferenceVectors(String[][] transforms, int pWidth, int pHeight, int w, int h, int xSkip, int ySkip, double max) throws InstantiationException, IllegalAccessException, NotEnoughDataPointsException, IllDefinedDataPointsException, ClassNotFoundException {
		ColorProcessor table = new ColorProcessor(
			(w + xSkip) * transforms.length - xSkip,
			(h + ySkip) * transforms.length - ySkip);
		ImagePlus impTable = new ImagePlus("Matrix", table);
		impTable.show();
		
		for (int i = 0; i < transforms.length; ++i) {
			CoordinateTransformList<mpicbg.models.CoordinateTransform> ct1 = createTransformList(i, transforms);
			for (int j = 0; j < transforms.length; ++j) {
				CoordinateTransformList<mpicbg.models.CoordinateTransform> ct2 = createTransformList(j, transforms);
		
				/* fit a simple linear model to compare with using some transferred samples */
				mpicbg.models.CoordinateTransform t = (mpicbg.models.CoordinateTransform)sampleModel2(ct2, ct1, invarianceModelClass, pWidth, pHeight);
				ct2.add(t);
		
				ImageProcessor ip = visualizeDifferenceVectors(
					w,
					h,
					pWidth,
					pHeight,
					ct1,
					ct2,
					max);
		
				table.copyBits(ip, (w + xSkip) * i, (h + ySkip) * j, Blitter.COPY);
				impTable.updateAndDraw();
			}
		}
		return impTable;
	}


	static ImagePlus showDifferenceVectorDistributions(String[][] transforms, int pWidth, int pHeight, int w, int h, int xSkip, int ySkip, double max) throws InstantiationException, IllegalAccessException, NotEnoughDataPointsException, IllDefinedDataPointsException, ClassNotFoundException {
		FloatProcessor table = new FloatProcessor(
			(w + xSkip) * transforms.length - xSkip,
			(h + ySkip) * transforms.length - ySkip);
		ImagePlus impTable = new ImagePlus("Matrix", table);
		impTable.show();
		
		
		for (int i = 0; i < transforms.length; ++i) {
			CoordinateTransformList<mpicbg.models.CoordinateTransform> ct1 = createTransformList(i, transforms);
			for (int j = 0; j < transforms.length; ++j) {
				CoordinateTransformList<mpicbg.models.CoordinateTransform> ct2 = createTransformList(j, transforms);
		
				/* fit a simple linear model to compare with using some transferred samples */
				mpicbg.models.CoordinateTransform t = (mpicbg.models.CoordinateTransform)sampleModel2(ct2, ct1, invarianceModelClass, pWidth, pHeight);
				ct2.add(t);
		
				ImageProcessor ip = visualizeDifferenceVectorDistribution(
					w,
					h,
					pWidth,
					pHeight,
					ct1,
					ct2,
					max);
		
				table.copyBits(ip, (w + xSkip) * i, (h + ySkip) * j, Blitter.COPY);
				impTable.updateAndDraw();
			}
		}
		return impTable;
	}

	static void drawCircles(ColorProcessor ip, String[][] transforms, int w, int h, int xSkip, int ySkip, double max) {
		for (int s = 1; s <= max; ++s) {
			ip.setColor(
				new Color(
					(float)(Math.min(1.0f, 2.0f * s / max)),
					(float)(Math.min(1.0f, 2.0f - 2 * s / max)),
					0.0f));
			for (int i = 0; i < transforms.length; ++i) {
				int x = (w + xSkip) * i;
				for (int j = i + 1; j < transforms.length; ++j) {
					int y = (h + ySkip) * j;
					ip.drawOval(
						(int)(w / 2 + x - Math.round(s * w / max / 2)),
						(int)(h / 2 + y - Math.round(s * h / max / 2)),
						(int)Math.round(s * w / max) + 1,
						(int)Math.round(s * h / max) + 1);
				}
			}
		}
	}

	static void drawLabels(ImagePlus imp, int offset, String[][] transforms, int w, int h, int xSkip, int ySkip) {
		ImageProcessor ip = imp.getProcessor().createProcessor(imp.getWidth() + offset, imp.getHeight() +offset);
		ip.copyBits(imp.getProcessor(), offset, offset, Blitter.COPY);
		ip.setColor(Color.WHITE);
		ip.setAntialiasedText(true);
		ip.setJustification(ImageProcessor.CENTER_JUSTIFY);
		Font font = new Font(Font.SANS_SERIF, Font.PLAIN, 14);
		ip.setFont(font);
		for (int i = 0; i < transforms.length; ++i) {
			int x = (w + xSkip) * i;
			ip.drawString(
				transforms[i][0],
				(int)(w / 2 + x + offset),
				20);
		}
		ip = ip.rotateRight();
//		ip.flipHorizontal();
		ip.setColor(Color.WHITE);
		ip.setAntialiasedText(true);
		ip.setJustification(ImageProcessor.CENTER_JUSTIFY);
		//ip.setFont(font.deriveFont(new AffineTransform(-1, 0, 0, 1, 0, 0)));
		ip.setFont(font);
		for (int i = 0; i < transforms.length; ++i) {
			int x = (w + xSkip) * i;
			ip.drawString(
				transforms[transforms.length - 1 - i][0],
				(int)(w / 2 + x),
				20);
		}
		ip = ip.rotateLeft();
//		ip.flipHorizontal();
		imp.setProcessor(ip);
	}


	
	public static HashMap<String, Object> exportTransform(String name, List<HashMap<String, String>> transform) {
		HashMap<String, Object> map = new HashMap<String, Object>();
		map.put("name", name);
		map.put("transform", transform);
		return map;
	}

	public static HashMap<String, String> exportTransform(CoordinateTransform ct) {
		HashMap<String, String> map = new HashMap<String, String>();
		String className = ct.getClass().getCanonicalName();
		if (className.equals("lenscorrection.NonLinearTransform"))
			className = "mpicbg.trakem2.transform.NonLinearCoordinateTransform";
		map.put("className", ct.getClass().getCanonicalName());
		map.put("dataString", ct.toDataString());
		return map;
	}

	public static ArrayList<PointMatch> samplePoints(Patch patch) {
		CoordinateTransform ct = patch.getFullCoordinateTransform();
		ArrayList<PointMatch> matches = new ArrayList<PointMatch>();
		double scaleX = (patch.getOWidth() - 1.0f) / 63.0f;
		double scaleY = (patch.getOHeight() - 1.0f) / 63.0f;
		for (int y = 0; y < 64; ++y) {
			double ys = scaleY * y;
			for (int x = 0; x < 64; ++x) {
				double xs = scaleX * x;
				Point p = new Point(new double[]{xs, ys});
				Point q = new Point(new double[]{xs, ys});
				ct.applyInPlace(p.getL());
				matches.add(new PointMatch(p, q));
			}
		}
		return matches;
	}
	
	public static List<String> findFiles(Path path, String[] fileExtensions) throws IOException {
		if (!Files.isDirectory(path)) {
			throw new IllegalArgumentException("path must be a directory.");
		}

		List<String> result;
		try (Stream<Path> flist = Files.list(path)) {
			result = flist.filter(p -> !Files.isDirectory(p))
					.map(p -> p.toString())
					.filter(f -> Arrays.stream(fileExtensions).anyMatch(f::endsWith)).collect(Collectors.toList());
		}
		return result;
	}
	
	public static ImagePlus ZMaxProjection(ImagePlus imp) {
		ZProjector zp = new ZProjector(imp);
		zp.setMethod(ZProjector.MAX_METHOD);
		if (imp.isHyperStack())
		{
			zp.setStopSlice(imp.getNSlices());
			zp.doHyperStackProjection(false);
		}
		else
			zp.doProjection();
		return zp.getProjection();
	}
	
	/**
	 * Interpolates between LUT values such that the result array has the
	 * specified number of colors.
	 * 
	 * @param baseLut
	 *            the LUT to interpolated
	 * @param nColors
	 *            the number of colors of the new LUT
	 * @return a nColors-by-3 array of color components
	 */
	public final static byte[][] interpolateLut(byte[][] baseLut, int nColors) {
		
		int n0 = baseLut.length;
		// allocate memory for new lut
		byte[][] lut = new byte[nColors][3];
		
		// linear interpolation of each color of new lut
		for (int i = 0; i < nColors; i++) {
			// compute color index in original lut
			float i0 = ((float) i) * n0 / nColors;
			int i1 = (int) Math.floor(i0);
			
			// the two surrounding colors
			byte[] col1 = baseLut[i1];
			byte[] col2 = baseLut[Math.min(i1 + 1, n0 - 1)];

			// the ratio between the two surrounding colors
			float f = i0 - i1;
			
			// linear interpolation of surrounding colors with cast
			lut[i][0] = (byte) ((1. - f) * (col1[0] & 0xFF) + f * (col2[0] & 0xFF));
			lut[i][1] = (byte) ((1. - f) * (col1[1] & 0xFF) + f * (col2[1] & 0xFF));
			lut[i][2] = (byte) ((1. - f) * (col1[2] & 0xFF) + f * (col2[2] & 0xFF));
		}
		
		return lut;
	}
	
	/**
	 * Creates a byte array representing the Fire LUT.
	 * 
	 * @param nColors
	 *            number of colors
	 * @return a nColors-by-3 array of color components.
	 */
	public final static byte[][] createFireLut(int nColors) {
		byte[][] lut = createFireLut();
		if (nColors != lut.length) 
			lut = interpolateLut(lut, nColors);
		return lut;
	}

	/**
	 * Creates a byte array representing the Fire LUT.
	 * 
	 * @return an array of color components.
	 */
	public final static byte[][] createFireLut() {
		// initial values
		int[] r = { 0, 0, 1, 25, 49, 73, 98, 122, 146, 162, 173, 184, 195, 207,
				217, 229, 240, 252, 255, 255, 255, 255, 255, 255, 255, 255,
				255, 255, 255, 255, 255, 255 };
		int[] g = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 14, 35, 57, 79, 101,
				117, 133, 147, 161, 175, 190, 205, 219, 234, 248, 255, 255,
				255, 255 };
		int[] b = { 0, 61, 96, 130, 165, 192, 220, 227, 210, 181, 151, 122, 93,
				64, 35, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 35, 98, 160, 223,
				255 };

		// create map
		byte[][] map = new byte[r.length][3];
		
		// cast elements
		for (int i = 0; i < r.length; i++) {
			map[i][0] = (byte) r[i];
			map[i][1] = (byte) g[i];
			map[i][2] = (byte) b[i];
		}
		
		return  map;
	}
	
	private static LUT createFireLUT(double maxVal)
	{
		byte[][] lut = createFireLut(256);
		byte[] red = new byte[256];
		byte[] green = new byte[256];
		byte[] blue = new byte[256];
		for (int i = 0; i < 256; i++)
		{
			red[i] 		= lut[i][0];
			green[i] 	= lut[i][1];
			blue[i] 	= lut[i][2];
		}
		IndexColorModel cm = new IndexColorModel(8, 256, red, green, blue);
		return new LUT(cm, 0, maxVal);
	}
	
	
	public static void main(String[] args)
	{
		Options options = new Options();

        Option in_op = new Option("i", "input", true, "input directory");
        in_op .setRequired(true);
        options.addOption(in_op);

        Option out_op = new Option("o", "output", true, "output directory");
        out_op.setRequired(true);
        options.addOption(out_op);

        Option param_op = new Option("p", "param", true, "path to parameter setting file");
        param_op.setRequired(true);
        options.addOption(param_op);
        
        Option scope_op = new Option("n", "scope", true, "scope name");
        scope_op.setRequired(true);
        options.addOption(scope_op);
        
        Option sample_op = new Option("s", "sample", true, "sample name");
        sample_op.setRequired(true);
        options.addOption(sample_op);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("LensDistCorrection", options);

            System.exit(1);
            return;
        }

        String dir_path = cmd.getOptionValue("input");
        String outdir = cmd.getOptionValue("output");
        String parampath = cmd.getOptionValue("param");
        String scope = cmd.getOptionValue("scope");
        String sample = cmd.getOptionValue("sample");
        String pname = scope + "_" + sample;
        
        System.out.println("input_dir: " + dir_path);
        System.out.println("output_dir: " + outdir);
        System.out.println("settings: " + parampath);
        System.out.println("project_name: " + pname);

        if (dir_path.endsWith(File.separator))
        	dir_path = dir_path.substring(0, dir_path.length() - 1);
        if (outdir.endsWith(File.separator))
        	outdir = outdir.substring(0, outdir.length() - 1);
        
		try
		{

			String jsontxt = "";
	        try
	        {
	        	jsontxt = new String ( Files.readAllBytes( Paths.get(parampath) ) );
	        } 
	        catch (IOException e) 
	        {
	            e.printStackTrace();
	        }
	        JSONObject jo = new JSONObject(jsontxt);
	        
	        int maxNumThreads = Runtime.getRuntime().availableProcessors();
	        
	        JSONObject montage_jo = jo.getJSONObject("montageLayers");
	        Align.ParamOptimize param = new Align.ParamOptimize();
			param.sift.initialSigma = (float)montage_jo.getDouble("initialSigma");
			param.sift.steps = montage_jo.getInt("steps");
			param.sift.minOctaveSize = montage_jo.getInt("minOctaveSize");
			param.sift.maxOctaveSize = montage_jo.getInt("maxOctaveSize");
			param.sift.fdSize = montage_jo.getInt("fdSize");
			param.sift.fdBins = montage_jo.getInt("fdBins");
			param.rod = (float)montage_jo.getDouble("rod");
			param.maxEpsilon = (float)montage_jo.getDouble("maxEpsilon");
			param.minInlierRatio = (float)montage_jo.getDouble("minInlierRatio");
			param.minNumInliers = montage_jo.getInt("minNumInliers");
			param.expectedModelIndex = montage_jo.getInt("expectedModelIndex");
			param.rejectIdentity = montage_jo.getBoolean("rejectIdentity");
			param.identityTolerance = (float)montage_jo.getDouble("identityTolerance");
			param.desiredModelIndex = montage_jo.getInt("desiredModelIndex");
			param.correspondenceWeight = (float)montage_jo.getDouble("correspondenceWeight");
			param.regularize = montage_jo.getBoolean("regularize");
			param.maxIterations = montage_jo.getInt("maxIterations");
			param.maxPlateauwidth = montage_jo.getInt("maxPlateauwidth");
			param.filterOutliers = montage_jo.getBoolean("filterOutliers");
			param.meanFactor = (float)montage_jo.getDouble("meanFactor");
	            
//			param.sift.initialSigma = 1.6f;
//			param.sift.steps = 3;
//			param.sift.minOctaveSize = 400;
//			param.sift.maxOctaveSize = 900;
//			param.sift.fdSize = 4;
//			param.sift.fdBins = 8;
//			param.rod = 0.92f;
//			param.maxEpsilon = 50.0f;
//			param.minInlierRatio = 0.0f;
//			param.minNumInliers = 20;
//			param.expectedModelIndex = 0;
//			param.rejectIdentity = false;
//			param.identityTolerance = 0.5f;
//			param.desiredModelIndex = 0;
//			param.correspondenceWeight = 1.0f;
//			param.regularize = false;
//			param.maxIterations = 2000;
//			param.maxPlateauwidth = 200;
//			param.filterOutliers = false;
//			param.meanFactor = 3.0f; 

			
			JSONObject align_jo = jo.getJSONObject("alignLayers");
			RegularizedAffineLayerAlignment.Param param2 = new RegularizedAffineLayerAlignment.Param(
					align_jo.getInt("SIFTfdBins"),//SIFTfdBins, 
					align_jo.getInt("SIFTfdSize"),//SIFTfdSize, 
					(float)align_jo.getDouble("SIFTinitialSigma"),//SIFTinitialSigma, 
					align_jo.getInt("SIFTmaxOctaveSize"),//SIFTmaxOctaveSize, 
					align_jo.getInt("SIFTminOctaveSize"),//SIFTminOctaveSize, 
					align_jo.getInt("SIFTsteps"),//SIFTsteps, 
					align_jo.getBoolean("clearCache"),//clearCache, 
					maxNumThreads,//maxNumThreadsSift,
					(float)align_jo.getDouble("rod"),//rod, 
					align_jo.getInt("desiredModelIndex"),//desiredModelIndex,
					align_jo.getInt("expectedModelIndex"),//expectedModelIndex, 
					(float)align_jo.getDouble("identityTolerance"),//identityTolerance,
					(float)align_jo.getDouble("lambda"),//lambda, 
					(float)align_jo.getDouble("maxEpsilon"),////maxEpsilon,
					align_jo.getInt("maxIterationsOptimize"),//maxIterationsOptimize,
					align_jo.getInt("maxNumFailures"),//maxNumFailures,
					align_jo.getInt("maxNumNeighbors"),//maxNumNeighbors, 
					maxNumThreads,//maxNumThreads, 
					align_jo.getInt("maxPlateauwidthOptimize"),//maxPlateauwidthOptimize,
					(float)align_jo.getDouble("minInlierRatio"),//minInlierRatio,
					align_jo.getInt("minNumInliers"),//minNumInliers,
					align_jo.getBoolean("multipleHypotheses"),//multipleHypotheses,
					align_jo.getBoolean("widestSetOnly"),//widestSetOnly,
					align_jo.getBoolean("regularize"),//regularize, 
					align_jo.getInt("regularizerIndex"),//regularizerIndex, 
					align_jo.getBoolean("rejectIdentity"),//rejectIdentity, 
					false//visualize
			);
			
//			RegularizedAffineLayerAlignment.Param param2 = new RegularizedAffineLayerAlignment.Param(
//					8,//SIFTfdBins, 
//					4,//SIFTfdSize, 
//					1.6f,//SIFTinitialSigma, 
//					1200,//SIFTmaxOctaveSize, 
//					400,//SIFTminOctaveSize, 
//					3,//SIFTsteps, 
//					true,//clearCache, 
//					maxNumThreads,//maxNumThreadsSift,
//					0.92f,//rod, 
//					0,//desiredModelIndex,
//					0,//expectedModelIndex, 
//					5.0f,//identityTolerance,
//					0.1f,//lambda, 
//					200.0f,////maxEpsilon,
//					1000,//maxIterationsOptimize,
//					5,//maxNumFailures,
//					5,//maxNumNeighbors, 
//					maxNumThreads,//maxNumThreads, 
//					200,//maxPlateauwidthOptimize,
//					0.0f,//minInlierRatio,
//					20,//minNumInliers,
//					true,//multipleHypotheses,
//					false,//widestSetOnly,
//					false,//regularize, 
//					1,//regularizerIndex, 
//					false,//rejectIdentity, 
//					false//visualize
//			);		
//			boolean propagateTransformBefore = false;
//			boolean propagateTransformAfter = false;
			

			JSONObject cd_jo = jo.getJSONObject("correctDistortion");
			CorrectDistortionFromSelectionParam p = new CorrectDistortionFromSelectionParam();
			p.sift.initialSigma = (float)cd_jo.getDouble("initialSigma");
			p.sift.steps = cd_jo.getInt("steps");
			p.sift.minOctaveSize = cd_jo.getInt("minOctaveSize");
			p.sift.maxOctaveSize = cd_jo.getInt("maxOctaveSize");
			p.sift.fdSize = cd_jo.getInt("fdSize");
			p.sift.fdBins = cd_jo.getInt("fdBins");
			p.rod = (float)cd_jo.getDouble("rod");
			p.maxNumThreadsSift = maxNumThreads;
			
			p.maxEpsilon = (float)cd_jo.getDouble("maxEpsilon");
			p.minInlierRatio = (float)cd_jo.getDouble("minInlierRatio");
			p.minNumInliers = cd_jo.getInt("minNumInliers");
			p.expectedModelIndex = cd_jo.getInt("expectedModelIndex");
			p.multipleHypotheses = cd_jo.getBoolean("multipleHypotheses");
			p.rejectIdentity = cd_jo.getBoolean("rejectIdentity");
			p.identityTolerance = (float)cd_jo.getDouble("identityTolerance");
			p.tilesAreInPlace = cd_jo.getBoolean("tilesAreInPlace");
			
			p.desiredModelIndex = cd_jo.getInt("desiredModelIndex");
			p.regularize = cd_jo.getBoolean("regularize");
			p.regularizerIndex = cd_jo.getInt("regularizerIndex");
			p.lambdaRegularize = (float)cd_jo.getDouble("lambdaRegularize");
			p.maxIterationsOptimize = cd_jo.getInt("maxIterationsOptimize");
			p.maxPlateauwidthOptimize = cd_jo.getInt("maxPlateauwidthOptimize");
			
			p.dimension = cd_jo.getInt("dimension");
			p.lambda = (float)cd_jo.getDouble("lambda");
			p.clearTransform = cd_jo.getBoolean("clearTransform");
			p.visualize = false;
			
//			CorrectDistortionFromSelectionParam p = new CorrectDistortionFromSelectionParam();
//			p.sift.initialSigma = 1.6f;
//			p.sift.steps = 3;
//			p.sift.minOctaveSize = 400;
//			p.sift.maxOctaveSize = 900;
//			p.sift.fdSize = 4;
//			p.sift.fdBins = 8;
//			p.rod = 0.92f;
//			p.maxNumThreadsSift = maxNumThreads;
//			
//			p.maxEpsilon = 5.0f;
//			p.minInlierRatio = 0.0f;
//			p.minNumInliers = 5;
//			p.expectedModelIndex = 1;
//			p.multipleHypotheses = true;
//			p.rejectIdentity = false;
//			p.identityTolerance = 5.0f;
//			p.tilesAreInPlace = true;
//			
//			p.desiredModelIndex = 0;
//			p.regularize = false;
//			p.regularizerIndex = 0;
//			p.lambdaRegularize = 0.01;
//			p.maxIterationsOptimize = 2000;
//			p.maxPlateauwidthOptimize = 200;
//			
//			p.dimension = 5;
//			p.lambda = 0.01f;
//			p.clearTransform = true;
//			p.visualize = false;
			
			JSONObject align2_jo = jo.getJSONObject("alignLayers2");
			RegularizedAffineLayerAlignment.Param param3 = new RegularizedAffineLayerAlignment.Param(
					align2_jo.getInt("SIFTfdBins"),//SIFTfdBins, 
					align2_jo.getInt("SIFTfdSize"),//SIFTfdSize, 
					(float)align2_jo.getDouble("SIFTinitialSigma"),//SIFTinitialSigma, 
					align2_jo.getInt("SIFTmaxOctaveSize"),//SIFTmaxOctaveSize, 
					align2_jo.getInt("SIFTminOctaveSize"),//SIFTminOctaveSize, 
					align2_jo.getInt("SIFTsteps"),//SIFTsteps, 
					align2_jo.getBoolean("clearCache"),//clearCache, 
					maxNumThreads,//maxNumThreadsSift,
					(float)align2_jo.getDouble("rod"),//rod, 
					align2_jo.getInt("desiredModelIndex"),//desiredModelIndex,
					align2_jo.getInt("expectedModelIndex"),//expectedModelIndex, 
					(float)align2_jo.getDouble("identityTolerance"),//identityTolerance,
					(float)align2_jo.getDouble("lambda"),//lambda, 
					(float)align2_jo.getDouble("maxEpsilon"),////maxEpsilon,
					align2_jo.getInt("maxIterationsOptimize"),//maxIterationsOptimize,
					align2_jo.getInt("maxNumFailures"),//maxNumFailures,
					align2_jo.getInt("maxNumNeighbors"),//maxNumNeighbors, 
					maxNumThreads,//maxNumThreads, 
					align2_jo.getInt("maxPlateauwidthOptimize"),//maxPlateauwidthOptimize,
					(float)align2_jo.getDouble("minInlierRatio"),//minInlierRatio,
					align2_jo.getInt("minNumInliers"),//minNumInliers,
					align2_jo.getBoolean("multipleHypotheses"),//multipleHypotheses,
					align2_jo.getBoolean("widestSetOnly"),//widestSetOnly,
					align2_jo.getBoolean("regularize"),//regularize, 
					align2_jo.getInt("regularizerIndex"),//regularizerIndex, 
					align2_jo.getBoolean("rejectIdentity"),//rejectIdentity, 
					false//visualize
			);
//			RegularizedAffineLayerAlignment.Param param3 = new RegularizedAffineLayerAlignment.Param(
//					8,//SIFTfdBins, 
//					4,//SIFTfdSize, 
//					1.6f,//SIFTinitialSigma, 
//					1200,//SIFTmaxOctaveSize, 
//					400,//SIFTminOctaveSize, 
//					3,//SIFTsteps, 
//					true,//clearCache, 
//					maxNumThreads,//maxNumThreadsSift,
//					0.92f,//rod, 
//					3,//desiredModelIndex,
//					0,//expectedModelIndex, 
//					5.0f,//identityTolerance,
//					0.01f,//lambda, 
//					200.0f,////maxEpsilon,
//					1000,//maxIterationsOptimize,
//					5,//maxNumFailures,
//					5,//maxNumNeighbors, 
//					maxNumThreads,//maxNumThreads, 
//					200,//maxPlateauwidthOptimize,
//					0.0f,//minInlierRatio,
//					20,//minNumInliers,
//					true,//multipleHypotheses,
//					false,//widestSetOnly,
//					true,//regularize, 
//					0,//regularizerIndex, 
//					false,//rejectIdentity, 
//					false//visualize
//			);
			
			JSONObject align3_jo = jo.getJSONObject("alignLayers3");
			RegularizedAffineLayerAlignment.Param param4 = new RegularizedAffineLayerAlignment.Param(
					align3_jo.getInt("SIFTfdBins"),//SIFTfdBins, 
					align3_jo.getInt("SIFTfdSize"),//SIFTfdSize, 
					(float)align3_jo.getDouble("SIFTinitialSigma"),//SIFTinitialSigma, 
					align3_jo.getInt("SIFTmaxOctaveSize"),//SIFTmaxOctaveSize, 
					align3_jo.getInt("SIFTminOctaveSize"),//SIFTminOctaveSize, 
					align3_jo.getInt("SIFTsteps"),//SIFTsteps, 
					align3_jo.getBoolean("clearCache"),//clearCache, 
					maxNumThreads,//maxNumThreadsSift,
					(float)align3_jo.getDouble("rod"),//rod, 
					align3_jo.getInt("desiredModelIndex"),//desiredModelIndex,
					align3_jo.getInt("expectedModelIndex"),//expectedModelIndex, 
					(float)align3_jo.getDouble("identityTolerance"),//identityTolerance,
					(float)align3_jo.getDouble("lambda"),//lambda, 
					(float)align3_jo.getDouble("maxEpsilon"),////maxEpsilon,
					align3_jo.getInt("maxIterationsOptimize"),//maxIterationsOptimize,
					align3_jo.getInt("maxNumFailures"),//maxNumFailures,
					align3_jo.getInt("maxNumNeighbors"),//maxNumNeighbors, 
					maxNumThreads,//maxNumThreads, 
					align3_jo.getInt("maxPlateauwidthOptimize"),//maxPlateauwidthOptimize,
					(float)align3_jo.getDouble("minInlierRatio"),//minInlierRatio,
					align3_jo.getInt("minNumInliers"),//minNumInliers,
					align3_jo.getBoolean("multipleHypotheses"),//multipleHypotheses,
					align3_jo.getBoolean("widestSetOnly"),//widestSetOnly,
					align3_jo.getBoolean("regularize"),//regularize, 
					align3_jo.getInt("regularizerIndex"),//regularizerIndex, 
					align3_jo.getBoolean("rejectIdentity"),//rejectIdentity, 
					false//visualize
			);
//			RegularizedAffineLayerAlignment.Param param4 = new RegularizedAffineLayerAlignment.Param(
//					8,//SIFTfdBins, 
//					4,//SIFTfdSize, 
//					1.6f,//SIFTinitialSigma, 
//					1200,//SIFTmaxOctaveSize, 
//					400,//SIFTminOctaveSize, 
//					3,//SIFTsteps, 
//					true,//clearCache, 
//					maxNumThreads,//maxNumThreadsSift,
//					0.92f,//rod, 
//					3,//desiredModelIndex,
//					0,//expectedModelIndex, 
//					5.0f,//identityTolerance,
//					0.01f,//lambda, 
//					50.0f,////maxEpsilon,
//					1000,//maxIterationsOptimize,
//					5,//maxNumFailures,
//					5,//maxNumNeighbors, 
//					maxNumThreads,//maxNumThreads, 
//					200,//maxPlateauwidthOptimize,
//					0.0f,//minInlierRatio,
//					20,//minNumInliers,
//					true,//multipleHypotheses,
//					false,//widestSetOnly,
//					true,//regularize, 
//					0,//regularizerIndex, 
//					false,//rejectIdentity, 
//					false//visualize
//			);
			
			String[] extensions = {"lsm", "LSM"};
			
			//new ImageJ();
			
			//read lsm files and generate mip images
			List<String> flist = findFiles(Paths.get(dir_path), extensions);		
			flist.sort(Comparator.naturalOrder());
			
			List<ImagePlus> mips = new ArrayList<ImagePlus>();
			for (int i = 0; i < flist.size(); i++) 
			{
				String path = flist.get(i);
				final ImagePlus[] impStack = BF.openImagePlus(path);

				System.out.println(
						"Number of images: " + impStack.length + " (this should always be 1), we ignore others");

				if (impStack.length > 1)
					throw new RuntimeException("More than one image was opened, please check the input carefully.");

				final ImagePlus imp = impStack[0];

				System.out.println("dimensions: " + imp.getStack().getProcessor(1).getWidth() + "x"
						+ imp.getStack().getProcessor(1).getHeight() + ", channels: " + imp.getNChannels()
						+ ", z-slices:" + imp.getNSlices() + ", timepoints: " + imp.getNFrames());

				imp.resetDisplayRange();
				
				ImagePlus mip_imp = ZMaxProjection(imp);
				HyperStackConverter.toStack(mip_imp);
				//ImagePlus mip_imp = ZProjector.run(imp, "max");
				//mip_imp.show();
				
				String fname = Paths.get(path).getFileName().toString();
				String prev_fname = "";
				if (i > 0)
					prev_fname = Paths.get(flist.get(i-1)).getFileName().toString();	
				
				if (i > 0 && fname.charAt(0) == prev_fname.charAt(0))
				{
					ImagePlus prev_mip = mips.get(mips.size()-1);
					ImageStack dst_stack = prev_mip.getStack();
					ImageStack src_stack = mip_imp.getStack();
					for (int j = 1; j <= src_stack.getSize(); j++)
						dst_stack.addSlice(src_stack.getProcessor(j).duplicate());
					prev_mip.setStack(dst_stack);
					mip_imp.close();
				}
				else
				{
					mips.add(mip_imp);
				}
				
				imp.close();
			}
			
			if (mips.size() == 0)
				throw new RuntimeException("mip creation failed.");
			
			int layernum = mips.get(0).getNSlices();
			int w = mips.get(0).getWidth();
			int h = mips.get(0).getHeight();
			List<ImageStack> layer_stacks = new ArrayList<ImageStack>();
			for (int i = 0; i < layernum; i++) 
				layer_stacks.add(new ImageStack(w, h));
			
			for (int i = 0; i < mips.size(); i++) 
			{
				ImagePlus mip = mips.get(i);
				ImageStack sstack = mip.getStack();
				for(int j = 0; j < layernum; j++)
				{
					layer_stacks.get(j).addSlice(sstack.getProcessor(j+1).duplicate());
				}
				mip.close();
			}
			
			List<ImagePlus> layers = new ArrayList<ImagePlus>();
			for (int i = 0; i < layernum; i++) 
			{
				ImagePlus layer_imp = new ImagePlus("layer"+i, layer_stacks.get(i));
				layers.add(layer_imp);
			}
			
			
			//save mip images
			//normalize local contrast brx 127 bry 127 stds 3.0 (all layers)
			String strage_dir = outdir+File.separator+pname;
			FileUtils.forceMkdir(new File(strage_dir));
			
			int brx = 127;
			int bry = 127;
			float stds = 3.0f;
			ArrayList<ArrayList<String>> layer_patch_paths = new ArrayList<ArrayList<String>>();
			for (int layer_id = 0; layer_id < layernum; layer_id++) 
			{
				ImagePlus layer = layers.get(layer_id);
				//layer.show();
				ImageStack sstack = layer.getStack();
				ArrayList<String> path_list = new ArrayList<String>();
				for (int i = 1; i <= sstack.getSize(); i++)
				{
					NormalizeLocalContrast.run(sstack.getProcessor(i), brx, bry, stds, true, true);
					String fname = String.format("layer_%02d_pos_%02d.tif", layer_id, i);
					ImagePlus tmp = new ImagePlus(fname, sstack.getProcessor(i).duplicate());
					FileSaver saver = new FileSaver(tmp);
					String fpath = strage_dir + File.separator + fname;
					saver.saveAsTiff(fpath);
					path_list.add(fpath);
				}
				layer.updateAndDraw();
				layer_patch_paths.add(path_list);
			}
			
			
			//create a new trakem project.
			ControlWindow.setGUIEnabled(false);
			
			Project project = Project.newFSProject("blank", null, strage_dir);
			LayerSet layerset = project.getRootLayerSet();
			for (int i = 0; i < layernum; i++)
				  layerset.getLayer(i, 1, true);
			project.getLayerTree().updateList(layerset);
			Display.updateLayerScroller(layerset);
			
			for (int i = 0; i < layernum; i++)
			{
				Layer layer = layerset.getLayer(i);
				ArrayList<String> path_list = layer_patch_paths.get(i);
				for (int s = 0; s < path_list.size(); s++)
				{
					Patch patch = Patch.createPatch(project, path_list.get(s));
					layer.add(patch);
				}
				layer.recreateBuckets();
			}
			
			//montage all layers. least square, translation.
			AlignTask.montageLayers(param, layerset.getLayers(), true, true, true, false, true);
			
			
			//Align layers. (least square)
			boolean propagateTransformBefore = false;
			boolean propagateTransformAfter = false;
			
			Rectangle box = null;
			HashSet< Layer > emptyLayers = new HashSet< Layer >();
			for ( final Iterator< Layer > it = layerset.getLayers().iterator(); it.hasNext(); )
			{
				/* remove empty layers */
				final Layer la = it.next();
				if ( !la.contains( Patch.class, true ) )
				{
					emptyLayers.add( la );
//					it.remove();
				}
				else
				{
					/* accumulate boxes */
					if ( null == box ) // The first layer:
						box = la.getMinimalBoundingBox( Patch.class, true );
					else
						box = box.union( la.getMinimalBoundingBox( Patch.class, true ) );
				}
			}
			
			new RegularizedAffineLayerAlignment().exec(param2, layerset.getLayers(), new HashSet<Layer>(), emptyLayers, box, propagateTransformBefore, propagateTransformAfter, null);
			
			
			//Auto resize canvas
			layerset.setMinimumDimensions();
			
			
			//Lens correction (All layers)
			for (int i = 0; i < layernum; i++)
			{
				p.firstLayerIndex = i;
				p.lastLayerIndex = i;
				final Layer layer = layerset.getLayer(i);
				ArrayList<Patch> patches = layer.getPatches(true);
				if (patches.size() > 0)
					DistortionCorrectionTask.run(p, patches, patches.get(0), layer);
			}
			
			
			//Align layers. least square	
			propagateTransformBefore = false;
			propagateTransformAfter = false;
			
			box = null;
			emptyLayers = new HashSet< Layer >();
			for ( final Iterator< Layer > it = layerset.getLayers().iterator(); it.hasNext(); )
			{
				/* remove empty layers */
				final Layer la = it.next();
				if ( !la.contains( Patch.class, true ) )
				{
					emptyLayers.add( la );
				}
				else
				{
					/* accumulate boxes */
					if ( null == box ) // The first layer:
						box = la.getMinimalBoundingBox( Patch.class, true );
					else
						box = box.union( la.getMinimalBoundingBox( Patch.class, true ) );
				}
			}
			
			new RegularizedAffineLayerAlignment().exec(param3, layerset.getLayers(), new HashSet<Layer>(), emptyLayers, box, propagateTransformBefore, propagateTransformAfter, null);
			
			
			//Align layers. least square		
			propagateTransformBefore = false;
			propagateTransformAfter = false;
			
			box = null;
			emptyLayers = new HashSet< Layer >();
			for ( final Iterator< Layer > it = layerset.getLayers().iterator(); it.hasNext(); )
			{
				/* remove empty layers */
				final Layer la = it.next();
				if ( !la.contains( Patch.class, true ) )
				{
					emptyLayers.add( la );
				}
				else
				{
					/* accumulate boxes */
					if ( null == box ) // The first layer:
						box = la.getMinimalBoundingBox( Patch.class, true );
					else
						box = box.union( la.getMinimalBoundingBox( Patch.class, true ) );
				}
			}
			
			new RegularizedAffineLayerAlignment().exec(param4, layerset.getLayers(), new HashSet<Layer>(), emptyLayers, box, propagateTransformBefore, propagateTransformAfter, null);
			
			
			//save trakem project
			project.saveAs(strage_dir + File.separator + pname + "_trakem_proj.xml", true);
			
			//String strage_dir = outdir+File.separator+pname;
			//Project project = Project.openFSProject(strage_dir + File.separator + pname + "_trakem_proj.xml");
			//LayerSet layerset = project.getRootLayerSet();
			
			//output coordinate transform
			Gson gson = new GsonBuilder().setPrettyPrinting().create();

			layerset.setMinimumDimensions();

			Rectangle topLeftBox = new Rectangle(100, 100);
			ArrayList<Patch> patches = new ArrayList<Patch>();

			for (Layer layer : layerset.getLayers()) {
				Collection<Displayable> displayables = layer.getDisplayables(Patch.class, topLeftBox);
				patches.add((Patch)displayables.iterator().next());
			}

			ArrayList<PointMatch> matches = new ArrayList<PointMatch>();

			for (Patch patch : patches)
				matches.addAll(samplePoints(patch));

			RigidModel2D model = new RigidModel2D();
			model.fit(matches);
			AffineModel2D affineModel = new AffineModel2D();
			affineModel.set(model.createAffine());

			ArrayList<HashMap<String, Object>> transformExports = new ArrayList<HashMap<String, Object>>();
			for (int i = 0; i < patches.size(); ++i) {
				Patch patch = patches.get(i);
				@SuppressWarnings("unchecked")
				final CoordinateTransformList< CoordinateTransform > ctl = (CoordinateTransformList< CoordinateTransform >) patch.getFullCoordinateTransform();
				List<CoordinateTransform> cts = ctl.getList(null);
				AffineModel2D affine = (AffineModel2D) cts.get(1);
				affine.preConcatenate(affineModel);
				List<HashMap<String, String>> maplist = new ArrayList<HashMap<String, String>>();
				maplist.add(exportTransform(ctl.get(0)));
				maplist.add(exportTransform(ctl.get(1)));
				HashMap<String, Object> export = exportTransform(
					String.format(format, scope, sample, lambdas[i]), maplist);
				transformExports.add(export);
			}
			
			String result_jsontxt = gson.toJson(transformExports);
			String jsonpath = outdir + File.separator + pname + ".json";
			Files.write(Paths.get(jsonpath), result_jsontxt.getBytes());
			System.out.println(result_jsontxt);
			
			
			//compare lenses
			ArrayList<ArrayList<String>> tr_lists = new ArrayList<ArrayList<String>>();
			ArrayList<String> temp_id_list = new ArrayList<String>();
			temp_id_list.add("Identity");
			temp_id_list.add("mpicbg.trakem2.transform.AffineModel2D");
			temp_id_list.add("1.0 0.0 0.0 1.0 0.0 0.0");
			tr_lists.add(temp_id_list);
			for (int i = 0; i < patches.size(); ++i) {
				Patch patch = patches.get(i);
				@SuppressWarnings("unchecked")
				final CoordinateTransformList< CoordinateTransform > ctl = (CoordinateTransformList< CoordinateTransform >) patch.getFullCoordinateTransform();
				List<CoordinateTransform> cts = ctl.getList(null);
				AffineModel2D affine = (AffineModel2D) cts.get(1);
				affine.preConcatenate(affineModel);
				ArrayList<String> templist = new ArrayList<String>();
				String label = String.format(format, scope, sample, lambdas[i]);
				String classname1 = ctl.get(0).getClass().getName();
				String ctstr1 = ctl.get(0).toDataString();
				String classname2 = ctl.get(1).getClass().getName();
				String ctstr2 = ctl.get(1).toDataString();
				templist.add(label);
				templist.add(classname1);
				templist.add(ctstr1);
				templist.add(classname2);
				templist.add(ctstr2);
				tr_lists.add(templist);
			}
			String[][] transforms = new String[tr_lists.size()][];
			String[] blankArray = new String[0];
			for(int i = 0; i < tr_lists.size(); i++) {
				transforms[i] = tr_lists.get(i).toArray(blankArray);
			}
			
			int iw = 256;
			int ih = 256;
			double max = 5;
			int pWidth = 1024;
			int pHeight = 1024;

			int ySkip = 4;
			int xSkip = 4;
			
			ImagePlus impVectors = showDifferenceVectors(transforms, pWidth, pHeight, iw, ih, xSkip, ySkip, max);
			ImagePlus impDists = showDifferenceVectorDistributions(transforms, pWidth, pHeight, iw, ih, xSkip, ySkip, max);

			impDists.setDisplayRange(0, 32);
			Thread.sleep(1000);
			impDists.setLut(createFireLUT(32.0));
			Thread.sleep(1000);
			new ImageConverter(impDists).convertToRGB();
			Thread.sleep(1000);
			impDists.getProcessor().snapshot();
			
//			{
//				FileSaver saver = new FileSaver(impDists);
//				saver.saveAsTiff(outdir + File.separator + pname + "_dists" + ".tif");
//				FileSaver saver2 = new FileSaver(impVectors);
//				saver2.saveAsTiff(outdir + File.separator + pname + "_vectors" + ".tif");
//			}
			
			ColorProcessor ip_src = (ColorProcessor)impVectors.getProcessor();
			ColorProcessor ip_dst = (ColorProcessor)impDists.getProcessor();
			
			for(int y = 0; y < ip_src.getHeight(); y++) {
				for(int x = y; x < ip_src.getWidth(); x++) {
					ip_dst.set(x, y, ip_src.get(x, y));
				}
			}
			
			drawCircles((ColorProcessor)impDists.getProcessor(), transforms, iw, ih, xSkip, ySkip, max);
			drawLabels(impDists, 26, transforms, iw, ih, xSkip, ySkip);
			
			FileSaver saver = new FileSaver(impDists);
			String compare_path = outdir + File.separator + pname + "_compare_lenses" + ".tif";
			saver.saveAsTiff(compare_path);
			

			System.out.println("Done");
			System.exit(0);
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit( 0 );
		}
	}
	
}
