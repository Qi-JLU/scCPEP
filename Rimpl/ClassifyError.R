# ClassifyError <- function(cellTypes_pred, cellTypes_test, cellTypes_train){

args <- commandArgs(trailingOnly = TRUE)
prd_type_dir <- args[1]
tst_type_dir <- args[2]
cellTypes_train_dir <- args[3]


#cell_annotations <- scan(annotation_dir, what=character())
cellTypes_pred <- readLines(prd_type_dir)
cellTypes_test <- readLines(tst_type_dir)
cellTypes_train <- readLines(cellTypes_train_dir)


errClass <- c("correct", "correctly unassigned",
              "intermediate", "incorrectly unassigned",
              "error assigned", "misclassified")

if (length(cellTypes_pred) != length(cellTypes_test)) {
    stop("wrong input")
}
train_ref <- unique(cellTypes_train)
res <- vapply(seq_len(length(cellTypes_pred)), function(i){
    if (cellTypes_test[i] %in% train_ref) {
        if (cellTypes_pred[i] %in% c("unassigned", "Unassigned")) {
            "incorrectly unassigned"
        } else if (cellTypes_pred[i] == "intermediate") {
            "intermediate"
        }else{
            if (cellTypes_test[i] == cellTypes_pred[i]) {
                "correct"
            }else if (grepl(cellTypes_test[i], cellTypes_pred[i])) {
                "intermediate"
            }
            else{
                "misclassified"
            }
        }
    }else{
        if (cellTypes_pred[i] %in% c("unassigned","Unassigned")) {
            "correctly unassigned"
        }else{
            "error assigned"
        }
    }
}, character(1L))

res <- factor(res, levels = errClass)
# print(as.matrix(res))
write(as.matrix(res), file.path(dirname(prd_type_dir), "res.txt"))
# return(res)
# }