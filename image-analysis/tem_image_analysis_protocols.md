# TEM image analysis protocols

Step-by-step procedures for the quantifications performed on transmission
electron microscopy (TEM) images in this study. Lipid droplets are measured at
lower magnification, where a single field contains enough droplets to quantify;
mitochondria and cristae are measured at higher magnification. These are
interactive workflows in Fiji/ImageJ and napari rather than scripts, and are
provided so that the measurements can be reproduced.

---

## 1. Lipid droplet area (TEM, 570x)

**Software:** Fiji/ImageJ v1.54
**Input:** TEM images (TIFF) acquired at 570x
**Scale:** 5 µm = 288.0156 pixels (57.6031 pixels/µm)

1. Open the image in Fiji and set the spatial scale (`Analyze > Set Scale`) to
   57.6031 pixels/unit before any selection is made, so that measurements are
   returned directly in µm².
2. Switch the selection tool to the Selection Brush Tool (right-click the oval
   selection icon in the toolbar). Left-click the icon to adjust the brush size.
3. For each droplet: clear any existing selection (`Ctrl+Shift+A`), outline the
   droplet with the brush, then press `T` to add it to the ROI Manager.
4. Repeat until every droplet in the field has been added.
5. In the ROI Manager, click `Deselect`, then `Measure`.
6. Export the Area column for downstream analysis.

**Notes.** Droplets are outlined manually; investigators performing the
outlining should be blinded to group identity, and the same magnification and
brush settings should be used across all images in a comparison. Record the
number of fields measured per animal and the number of animals per group, and
average per animal before statistical comparison.

---

## 2. Mitochondrial segmentation and cristae quantification (TEM, 1100x)

**Software:** napari with empanada-napari (MitoNet_v1_mini), Fiji/ImageJ with
Labkit, CLIJ2 and BioVoxxel 3D Box
**Input:** TEM images (TIFF) acquired at 1100x
**Scale:** 2 µm = 104 pixels (52 pixels/µm)

### 2.1 Mitochondrial instance and semantic masks (napari)

1. Activate the environment and launch napari:
   ```
   mamba activate empanada
   napari
   ```
2. Load the TIFF image.
3. `Plugins > empanada-napari > 2D Inference (Parameter Testing)`.
4. Select model `MitoNet_v1_mini` and run 2D inference. This produces the
   instance segmentation (referred to below as the *mask*).
5. Run 2D inference again with `Semantic only` selected. This produces the
   semantic segmentation (referred to below as the *label mask*).
6. Save both outputs.

### 2.2 Cristae segmentation (Fiji, Labkit)

1. Open the original TIFF in Fiji.
2. `Plugins > Labkit > Open Current Image With Labkit`.
3. Using the foreground brush, annotate representative cristae; using the
   background brush, annotate the remaining area. Corrections can be made with
   the eraser.
4. Train the classifier, inspect the result, and refine the annotations until
   the segmentation is satisfactory.
5. `Segmentation > Calculate Entire Segmentation Result`, then export the result
   (`File > Export Segmentation Result as Image`).

**Notes.** The classifier is trained interactively and its output depends on the
annotations provided. Train a single classifier and apply it unchanged to all
images in a comparison, and keep the annotator blinded to group identity.
Record the number of images used for training.

### 2.3 Cleaning the mitochondrial mask (Fiji, CLIJ2)

1. Open the mask in Fiji.
2. `CLIJ2-Assistant > Exclude Labels On Edges` to remove objects touching the
   image border.
3. On the result, `Label processing > Exclude Labels Outside Size Range (CLIJ2)`.
   Set the maximum size to the upper limit and the minimum size to the value
   used consistently across the dataset (10000 px in this study) to remove
   background and debris.
4. Save the cleaned label image.

**Notes.** The minimum size threshold must be fixed before analysis and applied
identically to every image in every group; report the value in the methods.

### 2.4 Restricting cristae to mitochondria (Fiji)

1. Open the label mask and the cristae segmentation.
2. Convert the label mask to a binary 0/1 image: `Ctrl+Shift+T` to threshold,
   `Apply`, convert to mask (0/255), then `Process > Math > Subtract` 254, and
   set the display range to 0/1 (`Ctrl+Shift+C`). Save.
3. `Process > Image Calculator`: multiply the binary label mask by the cristae
   segmentation. The product retains only cristae falling inside segmented
   mitochondria.
4. Optionally invert the lookup table (`Image > Lookup Tables > Invert LUT`) for
   display.

### 2.5 Measurement

1. Select the cleaned label image from step 2.3.
2. `Plugins > BioVoxxel 3D Box > Labels to 2D ROI Manager`.
3. Export the resulting measurements.

---

## Reporting

For both protocols, state in the manuscript the magnification, the pixel scale,
the number of fields per animal and animals per group, whether analysis was
blinded, and any fixed thresholds. Where several fields are measured per animal,
average within animal before statistical testing, or use a model that accounts
for the nesting.
