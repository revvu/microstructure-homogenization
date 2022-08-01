# Reevu Adakroy
# July 25, 2022
# arrange fibers in hexagonal orientation with chosen parameters and inputted volume fraction

from PIL import Image, ImageDraw
import math
import sys; args = sys.argv[1:]

# make input checking more robust
def validate_input():
    if len(args) != 1:
        print('Invalid input.')
        exit()
    return float(args[0])

def place_circles(n):
        lst = []
        start_height = (math.sqrt(3)/4)*n
        start_width = 0
        # first row
        for i in range(3):
            lst+= [(start_width+n/2+i*n, start_height)]
        # second row
        for i in range(4):
            lst+= [(start_width+i*n, start_height+math.sqrt(3)/2*n)]
        # third row
        for i in range(3):
            lst+= [(start_width+n/2+i*n, start_height+math.sqrt(3)*n)]
        # fourth row
        for i in range(4):
            lst+= [(start_width+i*n, start_height+3*math.sqrt(3)/2*n)]

        return lst

def main():
    volume_fraction = validate_input()
    fiber_count = 12
    height = 300
    width = 260
    radius = math.sqrt(height*width*volume_fraction/math.pi/fiber_count)
    n = width/3

    # store placement as a list of circle center coordinates
    circle_lst = place_circles(n)
    image = Image.new('RGBA', (width, height))

    draw = ImageDraw.Draw(image)
    draw.rectangle((0,0,width,height),fill=(54, 79, 107))
    for circle in circle_lst:
        draw.ellipse((circle[0] - radius, circle[1] - radius, circle[0] + radius, circle[1] + radius), fill=(63, 193, 201))
    image.show()
    image.save('hexagonal.png')

def execute(vol_frac, time_saved):
    volume_fraction = vol_frac
    fiber_count = 12
    height = 150
    width = 173
    # two fibers cut off
    radius = math.sqrt(height*width*volume_fraction/math.pi/9)
    n = width/3

    # store placement as a list of circle center coordinates
    circle_lst = place_circles(n)
    image = Image.new('RGBA', (width, height))

    draw = ImageDraw.Draw(image)
    draw.rectangle((0,0,width,height),fill=(54, 79, 107))
    for circle in circle_lst:
        draw.ellipse((circle[0] - radius, circle[1] - radius, circle[0] + radius, circle[1] + radius), fill=(63, 193, 201))

    filename = "generated_hexagonal" + time_saved + ".png"
    image.save('static/images/'+filename)

    with open('static/Generated Microstructures/RUC/generated_microstructure_'+time_saved+'.txt', 'w') as file:
        print(height,file=file)
        print(width,file=file)
        print(radius,file=file)
        for circle in circle_lst:
            print(circle[0],circle[1],file=file)


    return height, width


if __name__ == '__main__': main()