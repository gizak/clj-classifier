(ns gmproj2.core
  (:use clojure.pprint)
  (:gen-class))


(def dataset-dir "/Users/Gizak/Desktop/proj2_cross_field/")

(load "util")


(defn velocity-vec [p0 p1] ""
  [(- (first p1) (first p0)) (- (last p1) (last p0))])


(defn get-direction [vel] ""
  (let [h (cond
           (> (first vel) 0) "right"
           (< (first vel) 0) "left"
           :else "")
        v (cond
           (> (last vel) 0) "bottom"
           (< (last vel) 0) "top"
           :else "")]
    (keyword (string/join "-" [v h]))))


(def direction-weight {:top-left 1
                       :top-right 0
                       :top- 0.5
                       :bottom-left 0
                       :bottom-right 0
                       :bottom- 0.5
                       :-left 0.5
                       :-right 0.5
                       :- 0})


(defn feature-zigzag [traj] ""
  (loop [acc 0 remain (rest (seq traj))]
    (if (<= (count remain) 1) (/ (identity acc) (* 1 (count traj)))
        (recur (+ acc
                  ((get-direction (velocity-vec (first remain) (second remain)))
                   direction-weight))
               (rest remain)))))


(defn feature-curve-calc [traj mat] ""
  (loop [calc 0 remain (rest (seq traj))] ;; omit (0, 0)
    (if (empty? remain)  (/ (identity calc) 1)
        (recur (+ calc
                  (lookup mat :at (first remain)) 1)
               (rest remain)))))


(defn bunch-curve-calc [bunch mat] ""
  (map #(feature-curve-calc % mat) bunch))


(defn bunch-zigzag [bunch] ""
  (map #(feature-zigzag %) bunch))


(defn sum-zigzag [] ""
  (let [trajs (read-all-traj-bunch)
        bunch-done-zigzag (map
                         (fn [i] (bunch-zigzag (nth trajs i)))
                         (range 0 8))]
    (apply map + bunch-done-zigzag)))


(defn sum-curve-calc [] ""
  (let [maps (read-all-maps)
        trajs (read-all-traj-bunch)
        bunch-done-calc (map
                         (fn [i] (bunch-curve-calc (nth trajs i) (nth maps i)))
                         (range 0 8))]
    (apply map + bunch-done-calc)))


(defn avg-curve-calc [] ""
  (map #(float (/ % 8)) (sum-curve-calc)))


(defn avg-zigzag [] ""
  (map #(float (/ % 8)) (sum-zigzag)))


(def init-params {:classes [0 1 2]
                  :mus [[0.15 0.22 0.35]
                        [500 900 3500]]
                  :sigmas [[0.05 0.06 0.02]
                           [100 300 600]]
                  :class-dist [0.5 0.2 0.3]})


(defn class-posteriori [params x0 x1] ""
  (let [ndist (fn [i]
                (* (nth (:class-dist params) i)
                   (* (normal-dist (nth (first (:mus params)) i)
                                   (nth (first (:sigmas params)) i) x0)
                      (normal-dist (nth (second (:mus params)) i)
                                   (nth (second (:sigmas params)) i) x1))))
        rs (map #(ndist %) (:classes params))
        sum (apply + rs)]
    (map #(/ % sum) rs)))


(defn class-inference [params x0 x1] ""
  (let [ps (class-posteriori params x0 x1)]
    (nth (:classes params) (find-max-index ps))))


(defn compute-ess [params l0 l1] ""
  (loop [r0 (seq l0)
         r1 (seq l1)
         res []]
    (if (empty? r0) (identity res)
        (recur (rest r0) (rest r1)
               (conj res (class-posteriori params (first r0) (first r1)))))))


(defn compute-cprob-parti [ess] ""
  (apply map + ess))


(defn compute-cprob [ess] ""
  (let [part (compute-cprob-parti ess)
        sum (apply + part)]
    (map #(/ % sum) part)))


(defn compute-mu [ess l0 l1] ""
  (let [x (map #(identity [%1 %2]) l0 l1)
        pi (compute-cprob-parti ess)
        comp-ij (fn [i j x0 x1]
                  [(* (nth (nth ess i) j) x0)
                   (* (nth (nth ess i) j) x1)])
        map-ij (fn [j] (map #(comp-ij % j (first (nth x %)) (second (nth x %))) (range (count x))))
        sum-ij (fn [j]
                 (apply map + (map-ij j)))
        proc-j (fn [j] (map #(/ % (nth pi j)) (sum-ij j)))]
   (map #(proc-j %) (range (count pi)))))



(defn compute-sigma [ess l0 l1 mu] ""
  (let [x (map #(identity [%1 %2]) l0 l1)
        pi (compute-cprob-parti ess)
        comp-ij (fn [i j x0 x1]
                  [(* (nth (nth ess i) j) (math/expt (- (nth (first mu) j) x0) 2))
                   (* (nth (nth ess i) j) (math/expt (- (nth (second mu) j) x1) 2))])
        map-ij (fn [j] (map #(comp-ij % j (first (nth x %)) (second (nth x %))) (range (count x))))
        sum-ij (fn [j]
                 (apply map + (map-ij j)))
        proc-j (fn [j] (map #(/ % (nth pi j)) (sum-ij j)))]
   (map #(proc-j %) (range (count pi)))))


(defn compute-params [ess l0 l1 old] ""
  (let [crt (fn [arg] (vec (apply map #(identity [%1 %2 %3]) arg)))
        mu (crt (compute-mu ess l0 l1))
        sigma (crt (compute-sigma ess l0 l1 (:mus old)))
        cprob (compute-cprob ess)]
     (hash-map :classes [0 1 2]
               :mus mu
               :sigmas sigma
               :class-dist cprob)))


(defn expectation-maximization [l0 l1] ""
  (loop [params init-params
         ess nil]
    (let [new-ess (compute-ess params l0 l1)
              new-params (compute-params new-ess l0 l1 params)]
      (if (= params new-params) (identity params)
          (recur new-params new-ess)))))


(defn run-em [] ""
  (expectation-maximization (avg-zigzag) (avg-curve-calc)))


(defn dump-eng-den [path] ""
  (spit path (string/join "" (map #(format "%.3f\n" %) (avg-curve-calc)))))


(defn dump-zigzag [path] ""
  (spit path (string/join "" (map #(format "%.3f\n" %) (avg-zigzag)))))


(defn -main
  "I don't do a whole lot ... yet."
  [& args]
  (pprint (run-em)))
