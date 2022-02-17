package org.janelia.saalfeldlab.confocallens;

import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.HyperStackConverter;
import ij.plugin.ZProjector;
import ini.trakem2.Project;
import ini.trakem2.display.Display;
import ini.trakem2.display.Layer;
import ini.trakem2.display.LayerSet;
import ini.trakem2.display.Patch;
import loci.plugins.BF;
import mpicbg.ij.plugin.NormalizeLocalContrast;
import mpicbg.trakem2.align.Align;
import mpicbg.trakem2.align.AlignTask;
import mpicbg.trakem2.align.RegularizedAffineLayerAlignment;

public class Automation {
	
	public static List<String> findFiles(Path path, String[] fileExtensions) throws IOException {
		if (!Files.isDirectory(path)) {
			throw new IllegalArgumentException("path must be a directory.");
		}

		List<String> result;
		try (Stream<Path> flist = Files.list(path)) {
			result = flist.filter(p -> !Files.isDirectory(p))
					.map(p -> p.toString().toLowerCase())
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
	
	
	public static void main(String[] args)
	{
		
		String dir_path = "/Users/kawase/Scope9_20201119_40X_1024X1024";
		
		String[] extensions = {"lsm"};
		
		try
		{
			new ImageJ();
			
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
			String outdir = "/Users/kawase/lens_test";
			int brx = 127;
			int bry = 127;
			float stds = 3.0f;
			ArrayList<ArrayList<String>> layer_patch_paths = new ArrayList<ArrayList<String>>();
			for (int layer_id = 0; layer_id < layernum; layer_id++) 
			{
				ImagePlus layer = layers.get(layer_id);
				layer.show();
				ImageStack sstack = layer.getStack();
				ArrayList<String> path_list = new ArrayList<String>();
				for (int i = 1; i <= sstack.getSize(); i++)
				{
					NormalizeLocalContrast.run(sstack.getProcessor(i), brx, bry, stds, true, true);
					String fname = String.format("layer_%02d_pos_%02d", layer_id, i);
					ImagePlus tmp = new ImagePlus(fname, sstack.getProcessor(i).duplicate());
					FileSaver saver = new FileSaver(tmp);
					String fpath = outdir + File.separator + fname;
					saver.saveAsTiff(fpath);
					path_list.add(fpath);
				}
				layer.updateAndDraw();
				layer_patch_paths.add(path_list);
			}
			
			
			//create a new trakem project.
			Project project = Project.newFSProject("blank", null, "/Users/kawase/lens_test");
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
			Align.ParamOptimize param = new Align.ParamOptimize();
			param.sift.initialSigma = 1.6f;
			param.sift.steps = 3;
			param.sift.minOctaveSize = 400;
			param.sift.maxOctaveSize = 900;
			param.sift.fdSize = 4;
			param.sift.fdBins = 8;
			param.rod = 0.92f;
			param.maxEpsilon = 50.0f;
			param.minInlierRatio = 0.0f;
			param.minNumInliers = 20;
			param.expectedModelIndex = 0;
			param.rejectIdentity = false;
			param.identityTolerance = 0.5f;
			param.desiredModelIndex = 0;
			param.correspondenceWeight = 1.0f;
			param.regularize = false;
			param.maxIterations = 2000;
			param.maxPlateauwidth = 200;
			param.filterOutliers = false;
			param.meanFactor = 3.0f;
			
			AlignTask.montageLayers(param, layerset.getLayers(), true, true, true, false, true);
			
			
			//Align layers. (least square)
			int maxNumThreads = Runtime.getRuntime().availableProcessors();
			RegularizedAffineLayerAlignment.Param param2 = new RegularizedAffineLayerAlignment.Param(
					8,//SIFTfdBins, 
					4,//SIFTfdSize, 
					1.6f,//SIFTinitialSigma, 
					1200,//SIFTmaxOctaveSize, 
					400,//SIFTminOctaveSize, 
					3,//SIFTsteps, 
					true,//clearCache, 
					maxNumThreads,//maxNumThreadsSift,
					0.92f,//rod, 
					0,//desiredModelIndex,
					0,//expectedModelIndex, 
					5.0f,//identityTolerance,
					0.1f,//lambda, 
					200.0f,////maxEpsilon,
					1000,//maxIterationsOptimize,
					5,//maxNumFailures,
					5,//maxNumNeighbors, 
					maxNumThreads,//maxNumThreads, 
					200,//maxPlateauwidthOptimize,
					0.0f,//minInlierRatio,
					20,//minNumInliers,
					true,//multipleHypotheses,
					false,//widestSetOnly,
					false,//regularize, 
					1,//regularizerIndex, 
					false,//rejectIdentity, 
					false//visualize
			);
					
			boolean propagateTransformBefore = false;
			boolean propagateTransformAfter = false;
			
			Rectangle box = null;
			final HashSet< Layer > emptyLayers = new HashSet< Layer >();
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
			
			//Align layers. least square
			
			//output coordinate transform
			
			System.out.println("Done");
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit( 0 );
		}
	}
	
}
