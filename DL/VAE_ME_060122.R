### VAE first trial, using the example in the DL book
### Muhammad Elsadany
### This is June 1st, 2022



library(keras)
library(tensorflow)

### VAE encoder network
img_shape <- c(28,28,1)
batch_size <- 16
latent_dim <- 2L
input_img <- layer_input(shape=img_shape)

x <- input_img %>% 
  layer_conv_2d(filters = 32, kernel_size=3, padding="same", activation="relu") %>%
  layer_conv_2d(filters = 64, kernel_size=3, padding="same", activation="relu", strides=c(2,2))%>%
  layer_conv_2d(filters = 64, kernel_size=3, padding="same", activation="relu") %>%
  layer_conv_2d(filters = 64, kernel_size=3, padding="same", activation="relu")

shape_before_flattening <- k_int_shape(x)
x <- x %>% 
  layer_flatten() %>%
  layer_dense(units=32, activation="relu")

z_mean <- x %>%
  layer_dense(units= latent_dim)

z_log_var <- x %>% 
  layer_dense(units= latent_dim)


### latent-space sampling function
sampling <- function(args) {
  c(z_mean, z_log_var) %<-% args
  epsilon <- k_random_normal(shape = list(k_shape(z_mean) [1], latent_dim), mean=0, stddev=1)
  z_mean + k_exp(z_log_var) * epsilon
}

z <- list(z_mean, z_log_var) %>% 
  layer_lambda(sampling)

# decoder part 
decoder_input <- layer_input(k_int_shape(z) [-1])

x <- decoder_input %>%
  layer_dense(units = prod(as.integer(shape_before_flattening[-1])),
              activation = "relu") %>%
  layer_reshape(target_shape = shape_before_flattening[-1]) %>%
  layer_conv_2d_transpose(filters = 32, kernel_size = 3, padding = "same",
                          activation =  "relu", strides = c(2,2)) %>%
  layer_conv_2d(filters = 1, kernel_size = 3, padding = "same",
                activation = "sigmoid")

decoder <- keras_model(decoder_input, x)
z_decoded <- decoder(z)


# Custom layer used to compute the VAE loss
library(R6)
# this is checked and ready
CustomVariationalLayer <- R6Class("CustomVariationalLayer",
                                  inherit = KerasLayer,
                                  public = list(
                                    vae_loss = function(x, z_decoded) {
                                      x <- k_flatten(x)
                                      z_decoded <- k_flatten(z_decoded)
                                      xent_loss <- metric_binary_crossentropy(x, z_decoded)
                                      k1_loss <- -5e-4 * k_mean(
                                        1 + z_log_var - k_square(z_mean) - k_exp(z_log_var),
                                        axis = -1L
                                      )
                                      k_mean(xent_loss + k1_loss)
                                    },
                                    call = function(inputs, mask = NULL) {
                                      x < inputs[[1]]
                                      z_decoded <- inputs[[2]]
                                      loss <- self$vae_loss(x, z_decoded)
                                      self$add_loss(loss, inputs = inputs)
                                      x
                                    }
                                  )
)
layer_variational <- function(object) {
  create_layer(CustomVariationalLayer, object, list())
}
y <- list(input_img, z_decoded) %>%
  layer_variational()


