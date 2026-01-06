#################################
#     Functions in  DenMark     #
#################################



#-------- Function buildgridpp --------#

# build a function to discretize the MPP with grids:
# Input: matrix with locations, gene expression at each point 
# Output: two vectors (Y) and (M); 2 matrics of Y and M for the analysis 
buildgridpp<- function(grid_size=grid_size, 
                       dataset=dataset){
  #========input===================
  # dataset: a dataframe with locations (x, y) and 1 gene expression as the colnames 
  # grid_size: grid resolution at each of the x- and y-axis
  #================================
  # the min and max locations
  colnames(dataset)<- c('x', 'y', 'gene')
  loc<- as.data.frame(cbind(dataset$x, dataset$y))
  colnames(loc)<- c('x', 'y')
  # sum up the grids at each grid 
  loc_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
  marks_100<- matrix(NA, nrow=grid_size, ncol=grid_size)
  
  # initial values 
  init_x <- min(loc$x)
  init_y <- min(loc$y)
  
  # the delta: 
  step_x <- (max(loc$x)-min(loc$x))/grid_size
  step_y <- (max(loc$y)-min(loc$y))/grid_size
  
  # fill in the two matrix: 
  for (i in 1:grid_size){
    for (j in 1:grid_size){
      x_bound_1 <- init_x + (i-1) * step_x
      x_bound_2 <- init_x + i * step_x
      y_bound_1 <- init_y + (j-1) * step_y
      y_bound_2 <- init_y + j * step_y
      # index for each grid
      pts_index <-   intersect(which(loc$x >= x_bound_1 & 
                                       loc$x <= x_bound_2),
                               which(loc$y >= y_bound_1 & 
                                       loc$y <=  y_bound_2))
      # sum up the gene expression 
      marks_100[i,j]<- sum(dataset$gene[pts_index])
      loc_100[i,j] <- length(loc[pts_index,][,1])
    }
  }
  
  #======= return vectors ==== 
  # 1st element: 
  # the first column is the vectorized points 
  # the second column is the vectorized marks 
  # 2nd element: matrix for the points 
  # 3rd element: matrix for the marks 
  #===========================
  loc_mark_mat<- cbind(as.vector(loc_100), 
                       as.vector(marks_100))
  
  colnames(loc_mark_mat)<- c('point', 'mark')
  
  # summarize the results 
  loc_mark_res<- list(loc_mark_mat, 
                      loc_100, 
                      marks_100)
  
  return(loc_mark_res)
}



#---------- Function grid_pts_centroid ----------#

# Build another function for the locations of the centroids of those grids: 
# This is the function to sum up the characterisics we can get from the grids: 
# grid locations, grid_area, and the total number of grids 
grid_pts_centroid<- function(grid_size=grid_size, 
                             dataset=dataset){
  #========input===================
  # dataset: a dataframe with locations (x, y) and 1 gene expression as the colnames 
  # grid_size: grid resolution at each of the x- and y-axis
  #================================
  
  # the min and max locations
  loc<- as.data.frame(cbind(dataset$x, dataset$y))
  colnames(loc)<- c('x', 'y')
  # some key parameter for the bbox method
  max_coords_x <- max(loc$x)
  min_coords_x <- min(loc$x)
  max_coords_y <- max(loc$y)
  min_coords_y <- min(loc$y)
  
  # length of grid in two directions 
  len.grid.x <- ( max_coords_x - min_coords_x )/grid_size
  len.grid.y <- ( max_coords_y - min_coords_y )/grid_size
  
  # === The elements we care about ===  
  grid_area <- (len.grid.x * len.grid.y)
  
  # === The centroids of each grid === 
  grid_points <- expand.grid(x = seq(  min_coords_x-0.001,
                                       max_coords_x+0.001,
                                       length.out=grid_size), 
                             y = seq( min_coords_y-0.001,
                                      max_coords_y+0.001,
                                      length.out=grid_size))
  
  # === The number of total grid points 
  N <- nrow(grid_points)  # Total number of points
  
  
  #======= return list ==== 
  # the 1st element is the locations of the grids
  # the 2nd element is the grid area 
  # the 3rd element is the total number of the grid 
  #===========================
  grid_res_summary<- list(grid_points, 
                          grid_area, 
                          N)
  
  return(grid_res_summary)
}






