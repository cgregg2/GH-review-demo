from pyraf import iraf
from astropy.table import Table, Column
import numpy as np
import string, os, pyfits, glob, math, os.path #, pywcs
# TODO: update all to use astropy -- currently the pywcs import doesn't work

# pstamp: 
# updated 2006 Jan 19 to fix negative zero declination issue
def pstamp(images,coords,size,outdir='./',units='arcmin'):
    '''generates postage stamps from a list of images
      for a list of objects given in coords (format: ra dec name)
      with stampsize 'size' in units 'units'
      stamp images are rotated to have N-up, E-left.'''

    # carry size in degrees
    if units == 'arcsec':
        size = size/3600.0
    elif units == 'arcmin':
        size = size/60.0
    elif units == 'degrees':
        pass
    else:
        print 'Size must be either arcsec, arcmin, or degrees!'
        return

    # get the input
    imlist = parse_fname(images)
    coord_list = parse_coolist(coords)
    if len(imlist) == 0: 
        print 'Image list is empty; quitting'
        return
    if len(coord_list) == 0: 
        print 'Coordinate list is empty; quitting'
        return

    # loop through all the images
    for im in imlist:
      imbase = os.path.basename(im)
      rotangle, pxsc = get_wcs(im) 
      if pxsc == 0.0:
        print 'No WCS information in image %s' % im
        continue

      # loop through all the objects
      for ra, dec, objname in coord_list:
        outim = outdir + objname + '_'+ imbase
        if os.path.isfile(outim):
          print 'Output file named %s already exists' % outim
          continue

        # find the corners of the oversized stamp
        #  the factor of 0.71 allows for input image being rotated
        ra_big = ra + size*0.71
        ra_small = ra - size*0.71
        dec_big = dec + size*0.71
        dec_small = dec - size*0.71    
        x1, y1 = findpos(ra_big,dec_big,im)
        x2, y2 = findpos(ra_small,dec_small,im)
        x3, y3 = findpos(ra_big,dec_small,im)
        x4, y4 = findpos(ra_small,dec_big,im)
        xmax = int(round(max(x1,x2,x3,x4)))
        xmin = int(round(min(x1,x2,x3,x4)))
        ymax = int(round(max(y1,y2,y3,y4)))
        ymin = int(round(min(y1,y2,y3,y4)))

        # check to see if this is an HST image (need SCI extension)
        fimg = pyfits.open(im,'readonly')
        if len(fimg) > 1:
           extstr='[SCI]'
           hdr = fimg['SCI'].header
        else:
           hdr = fimg[0].header
           extstr=''

        # make sure none of the corners are outside the original image
        xsz = hdr['NAXIS1']
        ysz = hdr['NAXIS2']
        fimg.close()
        if (xmin,ymin)  < (0,0) or (xmax,ymax) > (xsz,ysz):
          print 'Coordinate %f %f falls outside image %s' %(ra,dec, im)
          continue


        # now make the oversized stamp and rotate it
        copystr = '[%d:%d,%d:%d]' % (xmin,xmax,ymin,ymax)
        if os.path.isfile('tmp.fits'):
          os.remove('tmp.fits')
        if os.path.isfile('tmpr.fits'):
          os.remove('tmpr.fits')
        iraf.imcopy(im+extstr+copystr,'tmp.fits',verbose='no')
        iraf.rotate('tmp.fits','tmpr.fits',rotangle,verbose='no')  # rotate to N-up

	# now trim to right size
        # here we let the object end up uncentered if necessary
        x, y = findpos(ra,dec,'tmpr.fits')  
        xsz, ysz = string.split(iraf.hselect('tmpr.fits','naxis[12]','yes',Stdout=1)[0],'\t')
        xsz, ysz = int(xsz), int(ysz)
        npix = size/pxsc
        xmin, xmax = x-npix*0.5, x+npix*0.5
        ymin, ymax = y-npix*0.5, y+npix*0.5
        if xmin < 1: xmin = 1
        if ymin < 1: ymin = 1
        if xmax > xsz: xmax = xsz
        if ymax > ysz: ymax = ysz

        copystr = '[%d:%d,%d:%d]' % (xmin,xmax,ymin,ymax)
        iraf.imcopy('tmpr.fits'+copystr,outim,verbose='no')
        print 'Extracted object %s from image %s' % (objname, im) # report success

#      os.remove('tmp.fits')
#      os.remove('tmpr.fits')

def parse_fname(name_list):
    """an IRAF-like input filename template parser"""
    if (len(name_list) > 0 and name_list[0] == '@'):
        try:                                                 # names are in a file
            fd = open(name_list[1:])                         #  one on each line
            text = fd.read()
            kw_list = text.split()
            if kw_list == []: kw_list = ['']
            fd.close()
        except IOError:
            print "Warning: file %s not found, assume null list" % name_list[1:]
            kw_list = ['']
    else:                                                                 
        kw_list = []
        _list = name_list.split(',')
	if len(_list) <= 0: kw_list = ['']
        for i in range(len(_list)):
            kw_list = kw_list + glob.glob(_list[i].strip())

    return kw_list


def parse_coolist(name_list):
    """parse the coordinate list for pstamp"""
    coo_list = []
    if (len(name_list) > 0 and name_list[0] == '@'):
        try:                                                 # objects are in a file
            fd = open(name_list[1:])                         #  one per line
            for line in fd.readlines():
              if line[0] == '#':
                 continue
              objinf = objtup(line)
              if len(objinf) == 3:
                 coo_list.append(objinf)
            fd.close()
        except IOError:
            pass
    else:      # just one object              
        objinf = objtup(name_list)
        if len(objinf) ==3:
           coo_list.append(objinf)
    return coo_list

# end of parse_coolist


def objtup(objinf):
    '''parse coordinate strings'''
     dat = objinf.split()
     if len(dat) != 3:
       return ()
     if string.find(dat[0],':') > 0:
       ra = string.split(dat[0],':') 
       dat[0] = 15.0*(float(ra[0])+float(ra[1])/60.0+float(ra[2])/3600.0)
       dec = string.split(dat[1],':') 
       if (float(dec[0])>0) or (dec[0][0]=='00') or (dec[0]=='+00'):
         dat[1] = float(dec[0])+float(dec[1])/60.0+float(dec[2])/3600.0
       else:
         dat[1] = float(dec[0])-float(dec[1])/60.0-float(dec[2])/3600.0
     else:
       dat[0] = float(dat[0])
       dat[1] = float(dat[1])         
     return (dat[0],dat[1],dat[2])

def findpos(ra,dec,im):
    '''locate the pixel coordinates corrdspnding to ra, dec in im'''
     sysstr = 'sky2xy %s %s %s' % (im, str(ra),str(dec))
     outstr = os.popen(sysstr).readlines()[0]
     x,y = string.split(outstr)[4],string.split(outstr)[5]
     if string.find(outstr,'off') > 0:
       x = -1
     return(float(x),float(y))

 def get_wcs(im):
     '''given an image, returns PA in degrees
         and pixel size (assumed to be square) in degrees'''
     fimg = pyfits.open(im,'readonly')
     if len(fimg) > 1:
        hdr = fimg['SCI'].header
     else:
        hdr = fimg[0].header
     if hdr.has_key('CROTA2') and hdr.has_key('CDELT2'): #  use CROTA1 and CDELT keywords if there
        cr2 = 0.0 - hdr['CROTA2']            
        pxsc = math.fabs(hdr['CDELT2'])
     elif hdr.has_key('CD1_2') and hdr.has_key('CD2_2'): #  or use CD matrix: lazy case assuming no skew
        cd12 = hdr['CD1_2']
        cd22 = hdr['CD2_2']
        cr2 = math.atan2(cd12,cd22)*180.0/math.pi  
        pxsc = math.sqrt(cd12**2+cd22**2)
     else:
        cr2 = 0.0
        pxsc = 0.0
     fimg.close()
     return(cr2,pxsc)

def sect_to_reg(infile,outfile):
    '''convert list of IRAF sections to ds9 region format'''
    fo = open(outfile,'w')
    fo.write('# Region file format: DS9 version 4.1\n')
    fo.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    fo.write('image\n')
    fi = open(infile,'r')
    for line in fi.readlines():
        xra, yra = string.split(line[1:-2],',')
        x1, x2 = float(string.split(xra,':')[0]), float(string.split(xra,':')[1])
        y1, y2 = float(string.split(yra,':')[0]), float(string.split(yra,':')[1])
        xcen = 0.5*(x1+x2)
        ycen = 0.5*(y1+y2)
        xw = x2-x1
        yw = y2-y1
        regtxt = 'box(%.1f,%.1f,%.0f,%.0f,0)\n' % (xcen, ycen, xw,yw)
        fo.write(regtxt)
    fi.close()
    fo.close()

def reg_to_sect(infile,outfile):
    '''turns a ds9 regions file into an IRAF-style section file (eg 1 3108 82 153)  
    only does rectangular sections, physical coordinates'''
    fo = open(outfile,'w')
    fi = open(infile,'r')
    for line in fi.readlines():
        if 'box' not in line:
            continue
        else:
            xcen, ycen, xw, yw, ang = string.split(line[4:-1],',')
#            if float(ang) != 0:
#                print "Can't deal with non-row/column regions, skipping %s" % line
#                continue
            x1 = float(xcen)-float(xw)/2.0
            x2 = float(xcen)+float(xw)/2.0
            y1 = float(ycen)-float(yw)/2.0
            y2 = float(ycen)+float(yw)/2.0
            regtxt = '%.1f %.1f %.1f %.1f\n' % (x1, x2, y1, y2)
        fo.write(regtxt)
    fi.close()
    fo.close()


def py_poly_phot(img_list, poly_list, output_sfx = '.ppy',coords='world'):
    '''do photmetry inside polygons'''
    #read list of polygon vertices 
    #same format as IRAF polyphot: one vertex per line, colon separates polygons
    polylist = []
    newpoly=[]
    fd = open(poly_list, 'r')
    for line in fd.readlines():
        if line[0] == ';':
            polylist.append(newpoly)
            newpoly=[]
        else:
            indat = string.split(line)
            x, y = float(indat[0]), float(indat[1])
            newpoly.append([x,y])
    fd.close()

    # loop through all the images
    for img in parse_fname(img_list):

	    #read image & header
	    im, hdr = pyfits.getdata(img, header=True)
	    # Parse the WCS keywords in the primary HDU
	    wcs = pywcs.WCS(hdr)
	    y, x = np.indices(im.shape).astype(float)
		
	    #loop over polygons
	    pindx = 0
	    tot_arr = np.zeros(len(polylist))
	    npix_arr = np.zeros(len(polylist))
	    mn_arr = np.zeros(len(polylist))
	    sig_arr = np.zeros(len(polylist))
	    unc_arr = np.zeros(len(polylist))
	    for poly in polylist:
                poly_img = convert_poly(poly, coords, wcs)
	        pindx +=1
                if len(poly) > 2:
		        msk = np.zeros(im.shape)     #create polygon mask
		        xmin, xmax, ymin, ymax = poly_corner(poly_img) # find the extent of the polygon
	                # decide how finely to divide up pixels based on extent
	                pix_area = (ymax-ymin)*(xmax-xmin)
	                if pix_area < 10:
	                    sd = 100
	                elif pix_area < 100:
	                    sd = 30
	                elif pix_area < 1000:
	                    sd = 20
	                else: 
	                    sd = 10
		        for x in range(xmin, xmax): # loop only over sub-region of image
		            for y in range(ymin, ymax):
		                msk[y,x] = pixel_inside_polygon(x,y,poly_img,subdiv=sd)
		
		        #compute masked statistics
		        ok_ind = np.where(np.logical_not(np.isnan(im)))
		        tot = (msk[ok_ind]*im[ok_ind]).sum()
		        npix = msk[ok_ind].sum() # note this was not quite correct before
		        mn = tot/npix
		        sig = math.sqrt(((msk[ok_ind]*im[ok_ind]*im[ok_ind]).sum()/npix) - mn*mn)
                        uncert_mn = sig/math.sqrt(npix)
                elif len(poly) == 1:
                    xpos = int(round(poly_img[0][0],0))
                    ypos = int(round(poly_img[0][1],0))
                    npix = 1
                    tot = im[ypos][xpos]
                    mn = im[ypos][xpos]
                    sig = 0.0
                    uncert_mn=0.0
                else:
                    print 'A line is not a polygon, silly'
                    continue
	        tot_arr[pindx-1]=tot
	        npix_arr[pindx-1]=npix
	        mn_arr[pindx-1]=mn
	        sig_arr[pindx-1]=sig
	        unc_arr[pindx-1]=uncert_mn
                print tot, npix,mn,sig,uncert_mn
                # end of loop over polygons

            # write the output to file
            imbase = os.path.basename(img)
            magout = imbase[:-5]+output_sfx
            if os.path.isfile(magout):
               os.unlink(magout)
            arrayascii.writeout(magout, [tot_arr, npix_arr, mn_arr, sig_arr, unc_arr]) # TODO: change to use np functions
            # end of loop over image files
    return


# grabbed online from http://www.ariel.com.au/a/python-point-int-poly.html
def point_inside_polygon(x,y,poly):
    '''determine if a point is inside a given polygon or not
    Polygon is a list of (x,y) pairs (the vertices)'''
    n = len(poly)
    inside =False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def pixel_inside_polygon(xpix, ypix, poly, subdiv=10):
    '''determine how much of a pixel centered at (xpix,ypix) is inside a polygon
    subdiv: pixel sub-division factor'''

    # recall that pixel centres are at integer values, so corners are half-integers
    # compute the coordinates of the pixel corners
    xmin = xpix-0.5
    xmax = xpix+0.5
    ymin = ypix-0.5
    ymax = ypix+0.5

    # check to see if any corners of the pixel are inside
    cnr1 = point_inside_polygon(xmin,ymin,poly)
    cnr2 = point_inside_polygon(xmin,ymax,poly)
    cnr3 = point_inside_polygon(xmax,ymin,poly)
    cnr4 = point_inside_polygon(xmax,ymax,poly)
    subsum = 0
    if (cnr1 and cnr2) and (cnr3 and cnr4): # if all 4 corners inside, fraction is 1
        subsum = subdiv*subdiv
    elif (cnr1 or cnr2) or (cnr3 or cnr4): # if at least one corner is inside, 
        xmin += 0.5/subdiv
        ymin += 0.5/subdiv
        for i in range(0, subdiv):  # divide  pixel up and loop
            for j in range(0, subdiv):
                subx = xmin + float(i)/subdiv
                suby = ymin + float(j)/subdiv
                subsum += point_inside_polygon(subx,suby,poly)

    # otherwise no corners inside polygon, fraction is zero
    infrac = float(subsum)/(subdiv*subdiv)
    return(infrac)
    

def poly_corner(poly): 
    '''take a list of (x,y) pairs and find the min,max values of both x and y'''
    # there's gotta be a more elegant way to do this, but hey
    x = np.zeros(len(poly))
    y = np.zeros(len(poly))
    for i in range(0,len(poly)):
        x[i] = poly[i][0]
        y[i] = poly[i][1]
    xmin =  int(math.floor(x.min())) 
    xmax = int(math.ceil(x.max()))
    ymin = int(math.floor(y.min())) 
    ymax = int(math.ceil(y.max()))
    return(xmin, xmax, ymin, ymax)

def convert_poly(poly, coords, wcs):
    '''convert the vertices of a polygon from one coordinate system to another'''
    newpoly=[]
    for vert in poly:
        if coords != 'world': #coords are in 1-indexed pixel coords
                x, y = vert[0]-1, vert[1]-1 # Convert  to 0-indexed pixel coordinates
        else: # Convert the sky coordinates to 0-indexed pixel coordinates
                pixcrd = wcs.wcs_sky2pix(np.array([vert]), 0) 
                x, y = pixcrd[0][0],pixcrd[0][1]
        newpoly.append([x,y])
    return(newpoly)

def parse_fname(name_list):
    """an IRAF-like input filename template parser"""
    if (len(name_list) > 0 and name_list[0] == '@'):
        try:                                                 # names are in a file
            fd = open(name_list[1:])                         #  one on each line
            text = fd.read()
            fd.close()
            kw_list = text.split()
            if kw_list == []: kw_list = ['']
            return kw_list
        except IOError:
            print "Warning: file %s not found, assume null list" % name_list[1:]
            return([''])
    else:                                                                 
        list = []
        _list = name_list.split(',')
        for i in range(len(_list)):
            list = list + glob.glob(_list[i].strip())
            
        return list

# NOT COMPLETE
def pymstat(img_glob, lower=None, upper=None, img_sect = None):
    """Pythonic version of IRAF imstat:
        produce image statistics for all science planes of an image,
        counting (but otherwise ignoring) NaNs

    Parameters
    ----------
    img_glob: string to glob for input images

    lower: float value below which to ignore pixel values (NOT IMPLEMENTED)

    upper: float value above which to ignore pixel values (NOT IMPLEMENTED)

    img_section: IRAF-style image-section over which to do stats (NOT IMPLEMENTED)
      (NB: uses 1-indexed pixel locations, [x_low:x_high, y_low: y_high])

    Outputs
    -------
    None (prints stats)"""


    imglist = glob.glob(img_glob)
    for img in imglist:
        hdulist = fits.open(img)
        img_planes = find_image_planes(hdulist)
        for plane in img_planes:
            dat = hdulist[plane].data
            valid = dat[~np.isnan(dat)]
            mean = valid.mean()
	    std = valid.std()
            med = valid.median()
            min=valid.min()
            max=valid.max()
            n_good = len(valid)
            n_nan = len(dat)-n_good
            # unfinished!!
            outstr = '{0} {:-f}; {:-f}'.format(plane_name,n_nan, n_good,mean,std,med,min,max) 
            print(outstr)
        hdulist.close()
    # end of loop over images
    return


def find_image_planes(hdulist):
    """
    Reads FITS hdulist to figure out which ones contain science data

    Parameters
    ----------
    hdulist: FITS hdulist

    Outputs
    -------
    img_plns: list of which indices in hdulist correspond to science data

    """
    n_hdu = len(hdulist)
    img_plns = []
    if n_hdu == 1: # if there is only one extension, then use that
        img_plns.append(0)
    else: # loop over all the extensions & try to find the right ones
        for extn in range(1,n_hdu):
            try: # look for 'EXTNAME' keyword, see if it's 'SCI'
                if 'SCI' in hdulist[extn].header['EXTNAME']:
                    img_plns.append(extn)
            except KeyError: # no 'EXTNAME', just assume we want this extension
                img_plns.append(extn)
    return(img_plns)

def transpose_table(tab_before, id_col_name='ID'):
    '''Returns a copy of tab_before (an astropy.Table) with rows and columns interchanged
        id_col_name: name for optional ID column corresponding to
        the column names of tab_before'''
    # contents of the first column of the old table provide column names for the new table
    # TBD: check for duplicates in new_colnames & resolve
    new_colnames=tuple(tab_before[tab_before.colnames[0]])
    # remaining columns of old table are row IDs for new table 
    new_rownames=tab_before.colnames[1:]
    # make a new, empty table
    tab_after=Table(names=new_colnames)
    # add the columns of the old table as rows of the new table
    for r in new_rownames:
        tab_after.add_row(tab_before[r])
    if id_col_name != '':
        # add the column headers of the old table as the id column of new table
        tab_after.add_column(Column(new_rownames, name=id_col_name),index=0)
    return(tab_after)
