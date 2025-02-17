use bullet_lib::{
    default::{inputs, loader, outputs, Loss, Trainer, TrainerBuilder},
    optimiser::{self, Optimiser},
    lr, wdl, Activation, LocalSettings, NetworkTrainer, TrainingSchedule, TrainingSteps,
};

type Output = outputs::MaterialCount<8>;
type Model = Trainer<optimiser::AdamWOptimiser, inputs::ChessBucketsMirrored, Output>;

const ID: &str = "approvers";
const DATA: &[&str] = &[r"..\dataset.bin"];

const THREADS: usize = 6;
const HIDDEN_SIZE: usize = 64;

const SCALE: f32 = 410.0;
const QA: i16 = 192;
const QB: i16 = 64;

const BATCHES: usize = 6850;
const SUPERBATCHES: usize = 350;

fn main() {
    stage1();
    stage2();
    stage3();
}

fn stage1() {
    let schedule = TrainingSchedule {
        net_id: format!("{ID}-stage1"),
        eval_scale: SCALE,
        steps: build_steps(),
        wdl_scheduler: wdl::ConstantWDL { value: 0.0 },
        lr_scheduler: lr::Warmup {
            inner: lr::LinearDecayLR {
                initial_lr: 0.001,
                final_lr: 0.001 * 0.3 * 0.3 * 0.3,
                final_superbatch: SUPERBATCHES,
            },
            warmup_batches: 800,
        },
        save_rate: SUPERBATCHES,
    };

    let mut trainer = build_trainer();

    trainer.set_optimiser_params(build_optimizer());

    let settings = build_settings();

    let data_loader = loader::DirectSequentialDataLoader::new(DATA);

    trainer.run(&schedule, &settings, &data_loader);
}

fn stage2() {
    let schedule = TrainingSchedule {
        net_id: format!("{ID}-stage2"),
        eval_scale: SCALE,
        steps: build_steps(),
        wdl_scheduler: wdl::ConstantWDL { value: 0.0 },
        lr_scheduler: lr::Warmup {
            inner: lr::CosineDecayLR {
                initial_lr: 0.001,
                final_lr: 0.001 * 0.3 * 0.3 * 0.3,
                final_superbatch: SUPERBATCHES,
            },
            warmup_batches: 200,
        },
        save_rate: SUPERBATCHES,
    };

    let mut trainer = build_trainer();

    trainer.optimiser_mut().load_weights_from_file(&checkpoint(1));

    trainer.set_optimiser_params(build_optimizer());

    let optimizer = optimiser::AdamWParams {
        min_weight: -1.0,
        max_weight:  1.0,
        ..build_optimizer()
    };

    trainer.optimiser_mut().set_params_for_weight("l1w", optimizer);
    trainer.optimiser_mut().set_params_for_weight("l1b", optimizer);

    let settings = build_settings();

    let data_loader = loader::DirectSequentialDataLoader::new(DATA);

    trainer.run(&schedule, &settings, &data_loader);
}

fn stage3() {
    let schedule = TrainingSchedule {
        net_id: format!("{ID}-stage3"),
        eval_scale: SCALE,
        steps: build_steps(),
        wdl_scheduler: wdl::ConstantWDL { value: 0.0 },
        lr_scheduler: lr::Warmup {
            inner: lr::LinearDecayLR {
                initial_lr: 0.001 * 0.3,
                final_lr: 0.001 * 0.3 * 0.3 * 0.3 * 0.3,
                final_superbatch: SUPERBATCHES,
            },
            warmup_batches: 200,
        },
        save_rate: SUPERBATCHES,
    };

    let mut trainer = build_trainer();

    trainer.optimiser_mut().load_weights_from_file(&checkpoint(2));

    trainer.set_optimiser_params(build_optimizer());

    let optimizer = optimiser::AdamWParams {
        min_weight: -f32::sqrt(2.0),
        max_weight:  f32::sqrt(2.0),
        ..build_optimizer()
    };

    trainer.optimiser_mut().set_params_for_weight("l1w", optimizer);
    trainer.optimiser_mut().set_params_for_weight("l1b", optimizer);

    let settings = build_settings();

    let data_loader = loader::DirectSequentialDataLoader::new(DATA);

    trainer.run(&schedule, &settings, &data_loader);
}

fn checkpoint(stage: usize) -> String {
    format!("checkpoints/{ID}-stage{stage}-{SUPERBATCHES}/optimiser_state/weights.bin")
}

fn build_steps() -> TrainingSteps {
    TrainingSteps {
        batch_size: 16384,
        batches_per_superbatch: BATCHES,
        start_superbatch: 1,
        end_superbatch: SUPERBATCHES,
    }
}

fn build_settings<'a>() -> LocalSettings<'a> {
    LocalSettings {
        test_set: None,
        threads: THREADS,
        output_directory: "checkpoints",
        batch_queue_size: 2048,
    }
}

fn build_optimizer() -> optimiser::AdamWParams {
    optimiser::AdamWParams {
        decay: 0.01,
        beta1: 0.9,
        beta2: 0.999,
        min_weight: f32::from(i8::MIN) / f32::from(QA),
        max_weight: f32::from(i8::MAX) / f32::from(QA),
    }
}

fn build_trainer() -> Model {
    TrainerBuilder::default()
        .quantisations(&[QA, QB])
        .optimiser(optimiser::AdamW)
        .loss_fn(Loss::SigmoidMSE)
        .input(inputs::ChessBucketsMirrored::new([0; 32]))
        .output_buckets(Output::default())
        .feature_transformer(HIDDEN_SIZE)
        .activate(Activation::SCReLU)
        .add_layer(1)
        .build()
}
