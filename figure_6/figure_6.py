from PIL import Image

#Read in each image
img_a = Image.open('./figure_6a.png')
img_b = Image.open('./figure_6b.png')
img_c = Image.open('./figure_6cd.png')

#Resample the scatter plot figure
img_a_i = img_a.resize(
    (int(float(img_a.size[0]*0.5)), int(float(img_a.size[1]*0.5))), Image.LANCZOS
)
img_c_i = img_c.resize(
    (int(float(img_c.size[0]*0.85)), int(float(img_c.size[1]*0.85))), Image.LANCZOS
)

#Make new image
out_size = (
    img_a_i.size[0] + img_b.size[0], img_a_i.size[1] + img_c_i.size[1] + 140 
)
output = Image.new("RGBA", out_size, color=(255, 255, 255))

#Paste images to output
output.paste(img_a_i, (0, 100))
output.paste(img_b, (img_a_i.size[0], 50))
output.paste(img_c_i, (100, int(img_a_i.size[1] + 140)))

#Save image
output.save('./figure_6.tiff')