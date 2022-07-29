# Reevu Adakroy
# July 28, 2022

# given the filename of a generated microstructure, appropriately generate an FVDAM input file
import math

def input_filename(time_saved):
    return 'static/Generated Microstructures/RUC/generated_microstructure_'+time_saved+'.txt'

def output_filename(time_saved):
    return 'static/Generated Microstructures/FVDAM Input/generated_microstructure_input_'+time_saved+'.fgm'

def parse_input(filename):
    with open(filename,'r') as file:
        lines = file.readlines()

    lines = [line.strip() for line in lines]
    height = int(lines[0])
    width = int(lines[1])
    radius = float(lines[2])

    circle_lst = [ tuple(float(a) for a in line.split()) for line in lines[3:]]

    return height, width, radius, circle_lst

def distance(x,y,circle_center): return math.sqrt((x-circle_center[0])**2+(y-circle_center[1])**2)


def build_matrix(height, width, radius, circle_lst):
    matrix = [['1' for _ in range(width)] for _ in range(height)]
    buffered_radius = 1.08*radius
    # for each circle traverse a diameter by diameter square around the circle
    for circle in circle_lst:
        for x in range(int(circle[0]-buffered_radius),int(circle[0]+buffered_radius)):
            if x<0: continue
            if x>=len(matrix[0]):break
            for y in range(int(circle[1]-buffered_radius),int(circle[1]+buffered_radius)):
                if y< 0: continue
                if y >= len(matrix): break
                if distance(x,y,circle) <= radius:
                    matrix[y][x]='2'

    return matrix

def main(time_saved):
    # read in file
    height, width, radius, circle_lst = parse_input(input_filename(time_saved))

    # discretize and generate material assignment matrix
    material_assignment = build_matrix(height, width, radius, circle_lst)

    # read through template and print matrix
    with open('fvdam_template.txt','r') as file:
        lines = file.readlines()

    with open(output_filename(time_saved), 'w') as file:
        line = 0
        while not "Nbeta	Ngama" in lines[line]:
            print(lines[line],file=file, end='')
            line+=1
        line+=1
        print("	"+str(height)+"	"+str(width)+"		Nbeta	Ngama", file=file)
        while line<len(lines):
            print(lines[line],file=file)
            line+=1

        for i in range(height):
            print(str(i+1)+'  0.000348346', file=file)
        print(file=file)

        for i in range(width):
            print('0.000348346',file=file)
        print(file=file)
        # print matrix
        for i in range(height, 0, -1):
            print(str(i)+'\t' + ' '.join(material_assignment[height-i]), file=file)
        print(file=file)

        return output_filename(time_saved)

if __name__ == '__main__': main('1659022426')
