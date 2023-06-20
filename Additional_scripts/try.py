
string1 = 'Hola:X-linked'
string2 = 'Hola:X-linked_recessive\x3bhola'
string3 = '.'

byteData = str.encode(string2)
string2_new = byteData.replace(b'\x3b', b':')
string2_decode = bytes.decode(string2_new)
print(string2_decode.split(':')[1])
print(string1.split(':')[1])
print(string3.split(':')[1])


