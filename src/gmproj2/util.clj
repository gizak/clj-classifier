(ns gmproj2.core
  (:require [clojure.string :as string]
            [clojure.contrib.math :as math]))

;; trim out the matrix
(defn read-matrix-from-file [^:String path]
  "Read the map/trajectories matrix from a given text file"
  (let [rows (-> path slurp (string/replace #"#.*?\n" "")
                 string/trim
                 string/split-lines)]
    (vec (map #(read-string (str "[" % "]")) rows))))


(defn lookup [matrix k [x y]]
  "Alias for matrix.At(x,y)"
  (nth (nth matrix (- x 1)) (- y 1)))


(defn split-trajectories [mat]
  "Find sub matrices from mat given delta row vector; Returns splited trajectories in a vector (i.e. 3-dimensions)"
  (loop [res [] remain (seq mat)]
    (if (empty? remain) (identity res)
      (if (= [0 0] (first remain))
        (recur (conj res [[0 0]]) (rest remain)) ;; assume start with [0 0]
        (recur (conj (pop res) (conj (last res) (first remain)))
               (rest remain))))))

;; trim out the matrix
(defn read-matrix-from-file [^:String path]
  "Read the map/trajectories matrix from a given text file"
  (let [rows (-> path slurp (string/replace #"#.*?\n" "")
                 string/trim
                 string/split-lines)]
    (vec (map #(read-string (str "[" % "]")) rows))))


(defn lookup [matrix k [x y]]
  "Alias for matrix.At(x,y)"
  (nth (nth matrix (- x 1)) (- y 1)))


(defn split-trajectories [mat]
  "Find sub matrices from mat given delta row vector; Returns splited trajectories in a vector (i.e. 3-dimensions)"
  (loop [res [] remain (seq mat)]
    (if (empty? remain) (identity res)
      (if (= [0 0] (first remain))
        (recur (conj res [[0 0]]) (rest remain)) ;; assume start with [0 0]
        (recur (conj (pop res) (conj (last res) (first remain)))
               (rest remain))))))


(defn read-all-mat [fname] ""
  (for [i (range 1 9)
    :let [path (str dataset-dir fname i ".txt")]]
    (read-matrix-from-file path)))


(defn read-all-maps [] ""
  (read-all-mat "field"))


(defn read-all-traj-bunch [] ""
  (map #(split-trajectories %) (read-all-mat "playOnField")))

(defn find-max-index [col] ""
  (loop [i 0 remain (seq col)]
    (if (empty? remain) (identity i)
        (if (> (first remain) (nth col i))
          (recur (- (count col) (count remain)) (rest remain))
          (recur i (rest remain))))))

(def e 3.14159265358979323846)
(def pi 2.7182818284590452353602874713526624977572)

(defn normal-dist [mu sigma x] ""
  (let [A (/ 1 (math/sqrt (* 2 3.14159265358979323846 sigma)))
        E (math/expt 2.7182818284590452353602874713526624977572
                     (/ (* -1 (math/expt (- x mu) 2))
                        (* 2 sigma)))]
    (* A E)))
