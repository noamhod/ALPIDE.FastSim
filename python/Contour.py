import cv2 
import numpy as np 

class Contour:
   '''
   see http://opencvpython.blogspot.com/2012/04/contour-features.html
   Provides detailed parameter informations about a contour
   Create a Contour instant as follows: c = Contour(pts,True)
   where pts is a simple list of 2D points (i.e. list of lists).

   Attributes:
   c.area -- gives the area of the region
   c.parameter -- gives the perimeter of the region
   c.moments -- gives all values of moments as a dict
   c.centroid -- gives the centroid of the region as a tuple (x,y)
   c.bounding_box -- gives the bounding box parameters as a tuple => (x,y,width,height)
   c.bx,c.by,c.bw,c.bh -- corresponds to (x,y,width,height) of the bounding box
   c.aspect_ratio -- aspect ratio is the ratio of width to height
   c.equi_diameter -- equivalent diameter of the circle with same as area as that of region
   c.extent -- extent = contour area/bounding box area
   c.convex_hull -- gives the convex hull of the region
   c.convex_area -- gives the area of the convex hull
   c.solidity -- solidity = contour area / convex hull area
   c.center -- gives the center of the ellipse
   c.majoraxis_length -- gives the length of major axis
   c.minoraxis_length -- gives the length of minor axis
   c.orientation -- gives the orientation of ellipse
   c.eccentricity -- gives the eccentricity of ellipse
   c.leftmost -- leftmost point of the contour
   c.rightmost -- rightmost point of the contour
   c.topmost -- topmost point of the contour
   c.bottommost -- bottommost point of the contour
   '''
   
   def __init__(self,lpts,doprint=False):
      self.lpts = lpts
      self.n    = len(lpts)
      if(self.n<5): self.Insert() ## fix lpts if len<4
      self.pts = np.array(self.lpts, np.float32)
      
      ## area
      self.area = cv2.contourArea(self.pts)
      if(doprint): print("area:",self.area)

      ## perimeter
      self.perimeter = cv2.arcLength(self.pts,True)
      if(doprint): print("perimeter:",self.perimeter)

      ## centroid
      self.moments = cv2.moments(self.pts)
      self.centroid = ()
      if self.moments['m00'] != 0.0:
         self.cx = self.moments['m10']/self.moments['m00']
         self.cy = self.moments['m01']/self.moments['m00']
         self.centroid = (self.cx,self.cy)
      else:
         self.centroid = "Region has zero area"
      if(doprint): print("centroid:",self.centroid)

      # bounding box
      self.bounding_box=cv2.boundingRect(self.pts)
      (self.bx,self.by,self.bw,self.bh) = self.bounding_box
      if(doprint): print("bounding_box:",self.bounding_box)

      # aspect ratio
      self.aspect_ratio = self.bw/float(self.bh)
      if(doprint): print("aspect_ratio:",self.aspect_ratio)

      # equivalent diameter
      self.equi_diameter = np.sqrt(4*self.area/np.pi)
      if(doprint): print("equi_diameter:",self.equi_diameter)

      # extent = contour area/boundingrect area
      self.extent = self.area/(self.bw*self.bh)
      if(doprint): print("extent:",self.extent)


      ### CONVEX HULL ###
      # convex hull
      self.convex_hull = cv2.convexHull(self.pts)
      if(doprint): print("convex_hull:",self.convex_hull)

      # convex hull area
      self.convex_area = cv2.contourArea(self.convex_hull)
      if(doprint): print("convex_area:",self.convex_area)

      # solidity = contour area / convex hull area
      self.solidity = self.area/float(self.convex_area)
      if(doprint): print("solidity:",self.solidity)


      if(self.n>=5):
         ### ELLIPSE  ###
         self.ellipse = cv2.fitEllipse(self.pts)
         if(doprint): print("ellipse:",self.ellipse)
   
         # center, axis_length and orientation of ellipse
         (self.center,self.axes,self.orientation) = self.ellipse
         if(doprint): print("ellipse(center,axes,orientation):",self.ellipse)
   
         # length of MAJOR and minor axis
         self.majoraxis_length = max(self.axes)
         self.minoraxis_length = min(self.axes)
         if(doprint): print("majoraxis_length:",self.majoraxis_length)
         if(doprint): print("minoraxis_length:",self.minoraxis_length)
   
         # eccentricity = sqrt( 1 - (ma/MA)^2) --- ma= minor axis --- MA= major axis
         self.eccentricity = np.sqrt(1-(self.minoraxis_length/self.majoraxis_length)**2)
         if(doprint): print("eccentricity:",self.eccentricity)


      ### CONTOUR APPROXIMATION ###
      self.contour_approx = cv2.approxPolyDP(self.pts,0.02*self.perimeter,True)
      if(doprint): print("contour_approx:",self.contour_approx)
      
      # ### EXTREME POINTS ###
      # # Finds the leftmost, rightmost, topmost and bottommost points
      # self.leftmost   = tuple(self.pts[self.pts[:,:,0].argmin()][0])
      # self.rightmost  = tuple(self.pts[self.pts[:,:,0].argmax()][0])
      # self.topmost    = tuple(self.pts[self.pts[:,:,1].argmin()][0])
      # self.bottommost = tuple(self.pts[self.pts[:,:,1].argmax()][0])
      # self.extreme    = (self.leftmost,self.rightmost,self.topmost,self.bottommost)
      # if(doprint): print("extreme(leftmost,rightmost,topmost,bottommost):",self.extreme)
      
      
   def Insert(self):
      if(self.n<5):
         pt0 = self.lpts[0]
         pt1 = self.lpts[1]
         ptnew = [(pt0[0]+pt1[0])/2, (pt0[1]+pt1[1])/2]
         # print("ptnew:",ptnew)
         self.lpts.insert(1,ptnew)
         self.n += 1
   

if __name__=='__main__':
   lpts = [[5.988352, 5.963775999999999], [6.986410666666666, 5.963775999999999], [6.986410666666666, 3.211264], [5.988352, 3.211264]]
   c = Contour(lpts,True)
   
   