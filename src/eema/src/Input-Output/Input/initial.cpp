#include "functions.h"

using namespace Eigen;

void initial_alloc() {

    int i, j, k ,l, m, n, iterator, ele_left, t_left, ret, covered;
    double result;

    nel_normal   = mesh[0].getNumElements();
    nnel_normal  = mesh[0].getNumNodesPerElement();
    nnode_normal = mesh[0].getNumNodes();
    sdof_normal  = nnode_normal * ndof;
    edof_normal  = nnel_normal * ndof;

    // Splitting of Threads

    unsigned threads_supported = std::thread::hardware_concurrency();

    number_of_threads = threads_supported;

    if (nel_normal <= number_of_threads) {
        start = new int[nel_normal];
        total = new int[nel_normal];

        for (iterator = 0; iterator < nel_normal; iterator++) {
            start[iterator] = iterator;
            total[iterator] = 1;
        }

        number_of_threads = nel_normal;
    }

    else {
        start = new int[number_of_threads];
        total = new int[number_of_threads];

        if (nel_normal % number_of_threads == 0) {
            for (iterator = 0; iterator < number_of_threads; iterator++) {
                start[iterator] = iterator * (nel_normal / number_of_threads);
                total[iterator] = nel_normal / number_of_threads;
            }
        }

        else {
            covered = 0;
            ele_left = nel_normal;
            t_left = number_of_threads;

            for (iterator = 0; iterator < number_of_threads - 1; iterator++) {
                result = ele_left/(double)(t_left);
                if (result + 0.5 > (int)result + 1)
                    ret = (int)result + 1;
                else
                    ret = (int)result;

                start[iterator] = covered;
                total[iterator] = ret;
                covered += ret;

                ele_left -= ret;
                t_left--;
            }

            start[iterator] = covered;
            total[iterator] = ele_left;
        }
    }


    // Allocating Memory

    I = MatrixXd::Identity(3, 3);
    points_normal  = guass_points(2);
    weights_normal = guass_weights(2);

    if (material_types_counter >= matTypeHigh) {
        matMap = new int[material_types_counter + 1];
        n = material_types_counter;
    }
    else {
        matMap = new int[matTypeHigh + 1];
        n = matTypeHigh;
    }

    for (iterator = 0; iterator < n; iterator++) {
        i = mat[iterator].getMatID();
        matMap[i] = iterator;
    }

    i = 0;
    n = 0;

    dndr_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndr_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndr_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndr_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndr_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    dnds_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dnds_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dnds_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dnds_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dnds_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    dndt_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndt_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndt_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndt_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndt_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    x_store = new double*[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        x_store[i] = new double[nnel_normal];
    }

    y_store = new double*[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        y_store[i] = new double[nnel_normal];
    }

    z_store = new double*[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        z_store[i] = new double[nnel_normal];
    }

    wtx_normal = new double[2];

    wty_normal = new double*[2];
    for (i = 0; i < 2; i++)
        wty_normal[i] = new double[2];

    wtz_normal = new double**[2];
    for (i = 0; i < 2; i++) {
        wtz_normal[i] = new double*[2];
        for (j = 0; j < 2; j++) {
            wtz_normal[i][j] = new double[2];
        }
    }

    jacobian_store = new double*****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        jacobian_store[i] = new double****[2];
        for (j = 0; j < 2; j++) {
            jacobian_store[i][j] = new double***[2];
            for (k = 0; k < 2; k++) {
                jacobian_store[i][j][k] = new double**[2];
                for (l = 0; l < 2; l++) {
                    jacobian_store[i][j][k][l] = new double*[ndof];
                    for (m = 0; m < ndof; m++) {
                        jacobian_store[i][j][k][l][m] = new double[ndof];
                    }
                }
            }
        }
    }

    invJacobian_store = new double*****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        invJacobian_store[i] = new double****[2];
        for (j = 0; j < 2; j++) {
            invJacobian_store[i][j] = new double***[2];
            for (k = 0; k < 2; k++) {
                invJacobian_store[i][j][k] = new double**[2];
                for (l = 0; l < 2; l++) {
                    invJacobian_store[i][j][k][l] = new double*[ndof];
                    for (m = 0; m < ndof; m++) {
                        invJacobian_store[i][j][k][l][m] = new double[ndof];
                    }
                }
            }
        }
    }

    det_store = new double***[nel_normal];
    for (j = 0; j < nel_normal; j++) {
        det_store[j] = new double**[2];
        for (k = 0; k < 2; k++) {
            det_store[j][k] = new double*[2];
            for (l = 0; l < 2; l++) {
                det_store[j][k][l] = new double[2];
            }
        }
    }

    dndx_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndx_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndx_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndx_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndx_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    dndy_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndy_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndy_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndy_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndy_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

    dndz_store = new double****[nel_normal];
    for (i = 0; i < nel_normal; i++) {
        dndz_store[i] = new double***[2];
        for (j = 0; j < 2; j++) {
            dndz_store[i][j] = new double**[2];
            for (k = 0; k < 2; k++) {
                dndz_store[i][j][k] = new double*[2];
                for (l = 0; l < 2; l++) {
                    dndz_store[i][j][k][l] = new double[nnel_normal];
                }
            }
        }
    }

  internalStressVariable1_prev_normal_store = new double*****[nel_normal];
  internalStressVariable2_prev_normal_store = new double*****[nel_normal];
  devInstantPK2Stress_prev_normal_store = new double*****[nel_normal];
  for (i = 0; i < nel_normal; i++) {
      internalStressVariable1_prev_normal_store[i] = new double****[2];
      internalStressVariable2_prev_normal_store[i] = new double****[2];
      devInstantPK2Stress_prev_normal_store[i] = new double****[2];
      for (j = 0; j < 2; j++) {
          internalStressVariable1_prev_normal_store[i][j] = new double***[2];
          internalStressVariable2_prev_normal_store[i][j] = new double***[2];
          devInstantPK2Stress_prev_normal_store[i][j] = new double***[2];
          for (k = 0; k < 2; k++) {
              internalStressVariable1_prev_normal_store[i][j][k] = new double**[2];
              internalStressVariable2_prev_normal_store[i][j][k] = new double**[2];
              devInstantPK2Stress_prev_normal_store[i][j][k] = new double**[2];
              for (l = 0; l < 2; l++) {
                  internalStressVariable1_prev_normal_store[i][j][k][l] = new double*[ndof];
                  internalStressVariable2_prev_normal_store[i][j][k][l] = new double*[ndof];
                  devInstantPK2Stress_prev_normal_store[i][j][k][l] = new double*[ndof];
                  for (m = 0; m < ndof; m++) {
                      internalStressVariable1_prev_normal_store[i][j][k][l][m] = new double[ndof];
                      internalStressVariable2_prev_normal_store[i][j][k][l][m] = new double[ndof];
                      devInstantPK2Stress_prev_normal_store[i][j][k][l][m] = new double[ndof];
                  }
              }
          }
      }
  }
  for (i = 0; i < nel_normal; i++) {
      for (j = 0; j < 2; j++) {
          for (k = 0; k < 2; k++) {
              for (l = 0; l < 2; l++) {
                  for (m = 0; m < ndof; m++) {
                    for (n = 0; n < ndof; n++) {
                        internalStressVariable1_prev_normal_store[i][j][k][l][m][n] = 0.0;
                        internalStressVariable2_prev_normal_store[i][j][k][l][m][n] = 0.0;
                        devInstantPK2Stress_prev_normal_store[i][j][k][l][m][n] = 0.0;
                    }
                  }
              }
          }
      }
  }

  internalStressVariable1_prev_centroid_store = new double**[nel_normal];
  internalStressVariable2_prev_centroid_store = new double**[nel_normal];
  devInstantPK2Stress_prev_centroid_store = new double**[nel_normal];
  for (i = 0; i < nel_normal; i++) {
      internalStressVariable1_prev_centroid_store[i] = new double*[ndof];
      internalStressVariable2_prev_centroid_store[i] = new double*[ndof];
      devInstantPK2Stress_prev_centroid_store[i] = new double*[ndof];
      for (j = 0; j < ndof; j++) {
          internalStressVariable1_prev_centroid_store[i][j] = new double[ndof];
          internalStressVariable2_prev_centroid_store[i][j] = new double[ndof];
          devInstantPK2Stress_prev_centroid_store[i][j] = new double[ndof];
      }
  }
  for (i = 0; i < nel_normal; i++) {
      for (j = 0; j < ndof; j++) {
        for (k = 0; k < ndof; k++) {
            internalStressVariable1_prev_centroid_store[i][j][k] = 0.0;
            internalStressVariable2_prev_centroid_store[i][j][k] = 0.0;
            devInstantPK2Stress_prev_centroid_store[i][j][k] = 0.0;
        }
      }
  }

}
