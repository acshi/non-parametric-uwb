#!/usr/bin/python3
import struct
import sys

def write_u8(f, val):
    f.write(struct.pack("B", val))

def write_u16(f, val):
    f.write(struct.pack("H", val))

def write_f32(f, val):
    f.write(struct.pack("<f", val))

def write_f64(f, val):
    f.write(struct.pack("<d", val))

for arg in sys.argv[1:]:
    with open(arg, "r") as read_f:
        with open("{}.binary".format(arg), "wb") as f:
            for line in read_f:
                parts = line.split()
                if len(parts) != 21:
                    break
                write_f64(f, float(parts[0]))
                write_u8(f, int(parts[1]))
                write_u8(f, int(parts[2]))
                write_f32(f, float(parts[3]))
                write_u16(f, int(parts[4]))
                write_u16(f, int(parts[5]))
                write_u16(f, int(parts[6]))
                write_f32(f, float(parts[7]))
                write_f32(f, float(parts[8]))
                write_f32(f, float(parts[9]))
                write_u16(f, int(parts[10]))
                write_u16(f, int(parts[11]))
                write_u16(f, int(parts[12]))
                write_f32(f, float(parts[13]))
                write_f32(f, float(parts[14]))
                write_f32(f, float(parts[15]))
                write_u16(f, int(parts[16]))
                write_u16(f, int(parts[17]))
                write_u16(f, int(parts[18]))
                write_f32(f, float(parts[19]))
                write_f32(f, float(parts[20]))
