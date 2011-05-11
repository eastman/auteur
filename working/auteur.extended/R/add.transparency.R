#generic plotting utility for adding transparency to a RGB color (specified by a 7 or 9 character string)
#author: R FITZ-JOHN 2010

add.transparency <-
function (col, alpha) {
    tmp <- col2rgb(col)/255
    rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = alpha)
}

